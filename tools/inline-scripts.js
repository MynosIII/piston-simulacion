// Copies the simulator scripts into Piston-Lab.html, so the page keeps
// working as a single offline file.
//   node tools/inline-scripts.js          update the HTML
//   node tools/inline-scripts.js --check  fail if the HTML is out of date
const fs = require('fs');
const path = require('path');
const root = path.join(__dirname, '..');
const htmlPath = path.join(root, 'Piston-Lab.html');
const files = ['Piston-Sim-Physics.js', 'Piston-Sim-UI.js'];

let html = fs.readFileSync(htmlPath, 'utf8');
const original = html;
for (const file of files) {
  const begin = `/* BEGIN ${file} */`, end = `/* END ${file} */`;
  const a = html.indexOf(begin), b = html.indexOf(end);
  if (a < 0 || b < a) throw new Error(`Missing ${begin} ... ${end} markers in Piston-Lab.html`);
  const source = fs.readFileSync(path.join(root, file), 'utf8').trimEnd();
  if (source.includes('</script')) throw new Error(`${file} cannot contain "</script"`);
  html = html.slice(0, a + begin.length) + '\n' + source + '\n' + html.slice(b);
}
if (process.argv.includes('--check')) {
  if (html !== original) { console.error('Piston-Lab.html is out of date: run node tools/inline-scripts.js'); process.exit(1); }
  console.log('Inline scripts are up to date.');
} else {
  fs.writeFileSync(htmlPath, html);
  console.log(html === original ? 'Already up to date.' : 'Piston-Lab.html updated.');
}
