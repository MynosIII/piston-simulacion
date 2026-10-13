# Pistón Lab

Simulación detallada de **un solo cilindro** de un motor de cuatro tiempos con encendido por chispa: el ciclo completo, el squish, el chorro de admisión, el swirl, las gotas de combustible, la película en las paredes, el frente de llama y el golpe (knock).

Es la hermana de [simulacion-motor](https://github.com/MynosIII/simulacion-motor) y sigue la misma filosofía: JavaScript puro, sin dependencias ni paso de compilación. Abrí `Piston-Lab.html` con doble clic y funciona sin conexión.

## Cómo se usa

- Arriba elegís un caso guardado: techo a dos aguas, cúpula, cámara descentrada, sin squish, inyección directa, carga estratificada, golpe, mezcla pobre, carga parcial, axisimétrico, etanol y GNC.
- La vista lateral muestra un campo en todo el corte, por defecto la composición. Con "Mitad derecha" podés poner otro campo en esa mitad para comparar. Los campos son: **composición** (mezcla azul, aire sin combustible azul oscuro, gases quemados grises, llama naranja, gotas blancas, película celeste), temperatura, vorticidad, velocidad, turbulencia, dosado o swirl.
- La **vista desde arriba** muestra el cilindro mirado por el eje: las dos válvulas de admisión, las dos de escape, la bujía y cada columna del corte promediada en la altura.
- La culata tiene sus **conductos de admisión y escape**, con el color del gas que contienen y rayas que avanzan con el caudal, y las válvulas se dibujan con su forma real (cara, asiento a 45°, vástago y guía), inclinadas según el techo.
- Los **trazadores** son marcadores que siguen al gas con su dispersión turbulenta: celestes los que entraron en este ciclo y naranjas el gas residual. Muestran cómo se mezcla la carga.
- Las flechas pueden mostrar la velocidad absoluta o la velocidad sin el movimiento que impone el pistón, que es donde se ven el tumble, el squish y los remolinos.
- A la derecha están la presión, el diagrama p–V, el calor liberado con la fracción quemada y el estado del combustible (vapor, líquido, película, quemado), más los valores del instante y los resultados del último ciclo.
- Abajo están todos los controles. Los de operación, inyección, válvulas y paredes cambian en vivo; los de geometría y combustible reinician la simulación.
- Espacio pausa o sigue, la flecha derecha avanza 1°.

La simulación corre en cámara lenta: un ciclo tarda unos segundos en pantalla según la velocidad elegida y la malla.

## Qué modela

**Geometría.** Hay dos cortes posibles:

- **Plano 2D** (por defecto): un corte que atraviesa el cilindro de lado a lado, con la admisión a la izquierda y el escape a la derecha. Ahí aparece el **tumble**. La profundidad del corte da el área real del pistón, así que masas y volúmenes son los del motor.
- **Axisimétrico**: un corte (r, y) alrededor del eje, con válvulas en anillo y **swirl**.

La culata puede ser de techo a dos aguas, cúpula, plana o con bowl en el pistón, y en el corte plano la cámara y la bujía se pueden descentrar. El área y el huelgo de squish son ajustables. La malla tiene filas en el recinto de la culata (cúpula), filas en el huelgo entre la cara del pistón y la culata, que se estiran con el pistón, y filas opcionales en un bowl del pistón. Las celdas fuera del perfil de la cámara quedan enmascaradas, así que las paredes son escalonadas. La relación de compresión es exacta.

**Flujo.** Formulación de bajo número de Mach: la presión termodinámica es uniforme y cambia en el tiempo, y una presión hidrodinámica empuja el gas entre celdas. En cada paso:

1. Fuentes: válvulas (orificio compresible con flujo inverso), inyección, gotas, película, chispa, llama y golpe.
2. Cada celda se comprime isentrópicamente a la nueva presión, elegida para que el gas llene exactamente el volumen nuevo.
3. Proyección: una ecuación de Poisson (Cholesky en banda) reparte entre vecinas el volumen que le sobra a cada celda. De ahí salen el squish, el chorro de admisión y la expansión de la llama.
4. Transporte conservativo con MUSCL y limitador minmod de masa, combustible, productos, energía, cantidad de movimiento, momento angular (swirl), turbulencia y fracción quemada.
5. Difusión implícita por líneas con viscosidad turbulenta, rozamiento y transferencia de calor en las paredes.

**Turbulencia.** Ecuación de k con producción por gradientes de velocidad, ley de pared, disipación con escala de mezcla limitada por las paredes y la altura de la cámara, y amplificación por compresión (distorsión rápida).

**Combustión.** En el corte plano la llama crece como un cilindro y no como una esfera, así que su velocidad lleva un factor fijo de 1.3 que devuelve duraciones de combustión de motor real. Un frente turbulento que avanza desde la bujía hacia el gas fresco a velocidad S_T = S_L (1 + 2 (u'/S_L)^0.7), con un factor de crecimiento del núcleo. La velocidad laminar S_L usa la correlación de Gülder con exponentes de Metghalchi-Keck para temperatura, presión y dilución por residuales. El frente se apaga en huelgos más finos que dos distancias de extinción, no entra en mezcla fuera del rango de inflamabilidad y deja combustible sin quemar en esos lugares. Los restos ricos que se encuentran con aire en gas caliente se oxidan después por mezcla.

**Combustible líquido.** Gotas lagrangianas con arrastre, calentamiento y evaporación (número de Spalding, Ranz-Marshall, Clausius-Clapeyron), acopladas en masa, energía y cantidad de movimiento al gas. Las que tocan una pared forman película, que se evapora hacia el gas o hierve si la pared está por encima de la temperatura de saturación. Con inyección en el puerto parte del combustible entra como gotas con la carga; con inyección directa el inyector central rocía un cono hueco en el ángulo elegido.

**Golpe (knock).** Integral de Livengood-Wu con el retardo de Douaud-Eyzat en la mezcla sin quemar. Cuando se completa, ese bolsón se quema de golpe. Se informa el ángulo, la masa que autoencendió y la frecuencia del primer modo acústico del cilindro (1.84 c / πB).

**Resultados por ciclo.** PMI neta y bruta, presión máxima, ángulos de 10, 50 y 90 % quemado, rendimiento indicado, eficiencia de combustión, calor a las paredes, llenado, λ atrapado, residuales, combustible sin quemar al abrir el escape, película máxima y golpe.

## Qué deja afuera

- La tercera dimensión: el corte plano tiene tumble pero no swirl, y el axisimétrico tiene swirl pero no tumble. Ninguno de los dos ve el flujo que cruza el corte.
- Las ondas de presión. El golpe se detecta y se estima su frecuencia, pero el pistoneo acústico no se simula.
- Química detallada, disociación, emisiones y hollín. Cada gas (aire, productos, vapor) tiene calor específico constante.
- La rotura primaria del chorro (las gotas nacen con su tamaño final), las salpicaduras al chocar con las paredes y los huelgos de los aros.

## Archivos

- `Piston-Lab.html`: la página, con los scripts incluidos para que funcione sola.
- `Piston-Sim-Physics.js`: el modelo físico. También se puede usar desde Node.
- `Piston-Sim-UI.js`: la vista, los gráficos y los controles.
- `tools/inline-scripts.js`: copia los dos scripts dentro del HTML. Corré `node tools/inline-scripts.js` después de editarlos (con `--check` solo verifica).
- `tests/physics.test.js`: pruebas de física. Corré `node tests/physics.test.js`.

Las pruebas verifican el solver en banda, que el corte plano genere tumble y queme como un motor, que el volumen respete la relación de compresión, la compresión isentrópica sin combustión, la primera ley en la parte cerrada del ciclo, las correlaciones de llama y retardo, la velocidad de squish contra la estimación de Heywood, la conservación del combustible entre vapor, gotas, película y quemado, valores razonables de un ciclo completo y que el golpe aparezca con avance excesivo.
