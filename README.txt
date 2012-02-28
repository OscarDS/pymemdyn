# Libreria ad-hoc para manejo de gromacs.

# Objetos de modelado. Estos tres archivos definen los objetos que se van a modelar.
- protein.py. Define los objetos Monomer, Ligand, Dimer y ProteinComplex. Estos objetos se inicializan con los archivos necesarios, y luego se pasan a otros objetos.
- membrane.py. Define la membrana.
- complex.py. Define el todo el complejo de proteina + membrana, y puede incluir cualquiera de los anteriores.

# Objetos auxiliares.
- queue.py. Sera el manejador de la cola al que se le pasaran los objetos de ejecucion. TODO.
- recipes.py. Lleva los pasos necesarios para ejecutar un modelado.
- utils.py. Lleva las funciones que los objetos anteriores utilizan a demanda. Por ejemplo, manipular archivos, copiar directorios, etc.

# Objeto de ejecucion.
- gromacs.py. Define los objetos Gromacs y Wrapper.
 * Gromacs carga los objetos de modelado, carga una receta de modelado y la ejecuta.
 * Wrapper es un proxy para los comandos gromacs. Se le pasa una entrada de la receta, y devuelve el comando a ejecutar.

# Archivo de ejemplo.
- example.py muestra como se utilizan las librerias anteriores.
