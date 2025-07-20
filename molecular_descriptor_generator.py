# ============================================================================
# ANÁLISIS DE DESCRIPTORES MOLECULARES 3D CON MÚLTIPLES MÉTODOS DE CONFORMACIÓN
# ============================================================================
# 
# Este script realiza un análisis exhaustivo de descriptores moleculares
# utilizando diferentes métodos de generación de conformaciones 3D y 
# calculando estadísticas de distribución para cada descriptor.

# ----------------------------------------------------------------------------
# IMPORTACIÓN DE LIBRERÍAS NECESARIAS
# ----------------------------------------------------------------------------

# Librerías principales para química computacional
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole

# Librerías para manipulación de datos y cálculos
import sys, py3Dmol
import pandas as pd
import numpy as np
import math

# Librería para cálculo de descriptores moleculares
from mordred import Calculator, descriptors

# Librería para visualización
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# CARGA Y PREPARACIÓN DE DATOS MOLECULARES
# ----------------------------------------------------------------------------

# Cargar dataset con estructuras moleculares desde archivo Excel
print("Cargando dataset molecular...")
df_original = pd.DataFrame(pd.read_excel('DATA.xlsx'))
lista_smiles = df_original["SMILES"]

print(f"Dataset cargado: {len(lista_smiles)} moléculas encontradas")

# Convertir notaciones SMILES a objetos moleculares de RDKit
moleculas_base = [Chem.MolFromSmiles(smile) for smile in lista_smiles]
moleculas_para_conformaciones = moleculas_base

# ----------------------------------------------------------------------------
# DEFINICIÓN DE FUNCIONES PARA GENERACIÓN DE CONFORMACIONES 3D
# ----------------------------------------------------------------------------

def generar_conformacion_uff(lista_moleculas):
    """
    Genera conformaciones 3D utilizando el campo de fuerza UFF (Universal Force Field)
    
    Args:
        lista_moleculas: Lista de objetos moleculares RDKit
        
    Returns:
        Lista de moléculas optimizadas con UFF
    """
    print("Generando conformaciones con campo de fuerza UFF...")
    moleculas_uff = []
    
    for molecula in lista_moleculas:
        # Añadir hidrógenos explícitos a la molécula
        mol_con_hidrogenos = Chem.AddHs(molecula)
        
        # Generar geometría 3D inicial sin conocimiento previo
        AllChem.EmbedMolecule(mol_con_hidrogenos, 
                            useBasicKnowledge=False, 
                            useExpTorsionAnglePrefs=False)
        
        # Verificar si UFF puede manejar todos los parámetros de la molécula
        if AllChem.UFFHasAllMoleculeParams(mol_con_hidrogenos):
            # Optimizar geometría con campo de fuerza UFF
            AllChem.UFFOptimizeMolecule(mol_con_hidrogenos)
            moleculas_uff.append(mol_con_hidrogenos)
    
    print(f"Conformaciones UFF generadas: {len(moleculas_uff)}")
    return moleculas_uff

def generar_conformacion_mmff(lista_moleculas):
    """
    Genera conformaciones 3D utilizando el campo de fuerza MMFF (Merck Molecular Force Field)
    
    Args:
        lista_moleculas: Lista de objetos moleculares RDKit
        
    Returns:
        Lista de moléculas optimizadas con MMFF
    """
    print("Generando conformaciones con campo de fuerza MMFF...")
    moleculas_mmff = []
    
    for molecula in lista_moleculas:
        # Añadir hidrógenos explícitos
        mol_con_hidrogenos = Chem.AddHs(molecula)
        
        # Generar geometría 3D inicial
        AllChem.EmbedMolecule(mol_con_hidrogenos, 
                            useBasicKnowledge=False, 
                            useExpTorsionAnglePrefs=False)
        
        # Verificar compatibilidad con MMFF
        if AllChem.MMFFHasAllMoleculeParams(mol_con_hidrogenos):
            # Optimizar con MMFF
            AllChem.MMFFOptimizeMolecule(mol_con_hidrogenos)
            moleculas_mmff.append(mol_con_hidrogenos)
    
    print(f"Conformaciones MMFF generadas: {len(moleculas_mmff)}")
    return moleculas_mmff

def generar_conformacion_etdg(lista_moleculas):
    """
    Genera conformaciones usando ETDG (Experimental-Torsion-angle preference with 
    Distance Geometry)
    
    Args:
        lista_moleculas: Lista de objetos moleculares RDKit
        
    Returns:
        Lista de moléculas con conformaciones ETDG
    """
    print("Generando conformaciones con método ETDG...")
    moleculas_etdg = []
    
    for molecula in lista_moleculas:
        mol_con_hidrogenos = Chem.AddHs(molecula)
        # Usar parámetros ETDG para geometría más realista
        AllChem.EmbedMolecule(mol_con_hidrogenos, AllChem.ETDG())
        moleculas_etdg.append(mol_con_hidrogenos)
    
    print(f"Conformaciones ETDG generadas: {len(moleculas_etdg)}")
    return moleculas_etdg

def generar_conformacion_kdg(lista_moleculas):
    """
    Genera conformaciones usando KDG (Knowledge-based Distance Geometry)
    
    Args:
        lista_moleculas: Lista de objetos moleculares RDKit
        
    Returns:
        Lista de moléculas con conformaciones KDG
    """
    print("Generando conformaciones con método KDG...")
    moleculas_kdg = []
    
    for molecula in lista_moleculas:
        mol_con_hidrogenos = Chem.AddHs(molecula)
        # Aplicar método KDG
        AllChem.EmbedMolecule(mol_con_hidrogenos, AllChem.KDG())
        moleculas_kdg.append(mol_con_hidrogenos)
    
    print(f"Conformaciones KDG generadas: {len(moleculas_kdg)}")
    return moleculas_kdg

def generar_conformacion_etkdg(lista_moleculas, version=1):
    """
    Genera conformaciones usando ETKDG (Experimental-Torsion-angle preference with 
    Knowledge-based Distance Geometry)
    
    Args:
        lista_moleculas: Lista de objetos moleculares RDKit
        version: Versión del algoritmo (1 o 2)
        
    Returns:
        Lista de moléculas con conformaciones ETKDG
    """
    print(f"Generando conformaciones con método ETKDGv{version}...")
    moleculas_etkdg = []
    
    for molecula in lista_moleculas:
        mol_con_hidrogenos = Chem.AddHs(molecula)
        
        # Seleccionar versión del algoritmo
        if version == 1:
            parametros = AllChem.ETKDG()
        elif version == 2:
            parametros = AllChem.ETKDGv2()
        else:
            print('Versión inválida. Use 1 o 2.')
            return None
            
        AllChem.EmbedMolecule(mol_con_hidrogenos, parametros)
        moleculas_etkdg.append(mol_con_hidrogenos)
    
    print(f"Conformaciones ETKDGv{version} generadas: {len(moleculas_etkdg)}")
    return moleculas_etkdg

# ----------------------------------------------------------------------------
# GENERACIÓN DE CONFORMACIONES CON DIFERENTES MÉTODOS
# ----------------------------------------------------------------------------

print("\n=== INICIANDO GENERACIÓN DE CONFORMACIONES ===")

# Generar conformaciones con cada método
conformaciones_uff = generar_conformacion_uff(moleculas_para_conformaciones)
conformaciones_etdg = generar_conformacion_etdg(moleculas_para_conformaciones)
conformaciones_etkdgv1 = generar_conformacion_etkdg(moleculas_para_conformaciones, version=1)
conformaciones_etkdgv2 = generar_conformacion_etkdg(moleculas_para_conformaciones, version=2)

# ----------------------------------------------------------------------------
# CONFIGURACIÓN DEL CALCULADOR DE DESCRIPTORES
# ----------------------------------------------------------------------------

# Definir número de iteraciones para análisis estadístico
numero_iteraciones = 50
lista_iteraciones = list(range(numero_iteraciones))
indices_moleculas = list(range(len(conformaciones_uff)))

print(f"\nConfiguración del análisis:")
print(f"- Número de iteraciones: {numero_iteraciones}")
print(f"- Número de moléculas: {len(conformaciones_uff)}")

# Crear calculador temporal para filtrar descriptores
calculador_temporal = Calculator(descriptors, ignore_3D=False)

# Lista específica de descriptores 3D de interés
nombres_descriptores_3d = [
    # Descriptores de superficie molecular polar/no polar
    "PNSA1","PNSA2","PNSA3","PNSA4","PNSA5",  # Superficie polar negativa
    "PPSA1","PPSA2","PPSA3","PPSA4","PPSA5",  # Superficie polar positiva
    "DPSA1","DPSA2","DPSA3","DPSA4","DPSA5",  # Diferencia de superficie polar
    "FNSA1","FNSA2","FNSA3","FNSA4","FNSA5",  # Superficie negativa fraccionada
    "FPSA1","FPSA2","FPSA3","FPSA4","FPSA5",  # Superficie positiva fraccionada
    "WNSA1","WNSA2","WNSA3","WNSA4","WNSA5",  # Superficie negativa ponderada
    "WPSA1","WPSA2","WPSA3","WPSA4","WPSA5",  # Superficie positiva ponderada
    
    # Descriptores topológicos de superficie
    "RNCS","RPCS","TASA","TPSA","RASA","RPSA",
    
    # Descriptores geométricos
    "GeomDiameter","GeomRadius","GeomShapeIndex","GeomPetitjeanIndex",
    
    # Descriptores gravitacionales
    "GRAV","GRAVH","GRAVp","GRAVHp",
    
    # Descriptores 3D-MoRSE (diversos tipos)
    "Mor01","Mor02","Mor03","Mor04","Mor05","Mor06","Mor07","Mor08","Mor09","Mor10",
    "Mor11","Mor12","Mor13","Mor14","Mor15","Mor16","Mor17","Mor18","Mor19","Mor20",
    "Mor21","Mor22","Mor23","Mor24","Mor25","Mor26","Mor27","Mor28","Mor29","Mor30",
    "Mor31","Mor32",
    
    # 3D-MoRSE ponderados por masa
    "Mor01m","Mor02m","Mor03m","Mor04m","Mor05m","Mor06m","Mor07m","Mor08m","Mor09m","Mor10m",
    "Mor11m","Mor12m","Mor13m","Mor14m","Mor15m","Mor16m","Mor17m","Mor18m","Mor19m","Mor20m",
    "Mor21m","Mor22m","Mor23m","Mor24m","Mor25m","Mor26m","Mor27m","Mor28m","Mor29m","Mor30m",
    "Mor31m","Mor32m",
    
    # 3D-MoRSE ponderados por volumen de van der Waals
    "Mor01v","Mor02v","Mor03v","Mor04v","Mor05v","Mor06v","Mor07v","Mor08v","Mor09v","Mor10v",
    "Mor11v","Mor12v","Mor13v","Mor14v","Mor15v","Mor16v","Mor17v","Mor18v","Mor19v","Mor20v",
    "Mor21v","Mor22v","Mor23v","Mor24v","Mor25v","Mor26v","Mor27v","Mor28v","Mor29v","Mor30v",
    "Mor31v","Mor32v",
    
    # 3D-MoRSE ponderados por electronegatividad de Sanderson
    "Mor01se","Mor02se","Mor03se","Mor04se","Mor05se","Mor06se","Mor07se","Mor08se","Mor09se","Mor10se",
    "Mor11se","Mor12se","Mor13se","Mor14se","Mor15se","Mor16se","Mor17se","Mor18se","Mor19se","Mor20se",
    "Mor21se","Mor22se","Mor23se","Mor24se","Mor25se","Mor26se","Mor27se","Mor28se","Mor29se","Mor30se",
    "Mor31se","Mor32se",
    
    # 3D-MoRSE ponderados por polarizabilidad
    "Mor01p","Mor02p","Mor03p","Mor04p","Mor05p","Mor06p","Mor07p","Mor08p","Mor09p","Mor10p",
    "Mor11p","Mor12p","Mor13p","Mor14p","Mor15p","Mor16p","Mor17p","Mor18p","Mor19p","Mor20p",
    "Mor21p","Mor22p","Mor23p","Mor24p","Mor25p","Mor26p","Mor27p","Mor28p","Mor29p","Mor30p",
    "Mor31p","Mor32p",
    
    # Momentos de inercia
    "MOMI-X","MOMI-Y","MOMI-Z","PBF"
]

print(f"- Descriptores 3D seleccionados: {len(nombres_descriptores_3d)}")

# Filtrar descriptores disponibles
descriptores_filtrados = []
for i, descriptor in enumerate(calculador_temporal.descriptors):
    if descriptor.__str__() in nombres_descriptores_3d:
       descriptores_filtrados.append(descriptor)

# Crear calculador final con descriptores filtrados
calculador_3d = Calculator(descriptores_filtrados, ignore_3D=False)

print(f"- Descriptores disponibles para cálculo: {len(descriptores_filtrados)}")

# ----------------------------------------------------------------------------
# CÁLCULO DE DESCRIPTORES PARA CADA ITERACIÓN Y MÉTODO
# ----------------------------------------------------------------------------

print("\n=== INICIANDO CÁLCULO DE DESCRIPTORES ===")

# Diccionario para almacenar resultados de cada iteración
resultados_descriptores = {}

# Calcular descriptores para cada iteración
for iteracion in lista_iteraciones:
    print(f"Procesando iteración {iteracion + 1}/{numero_iteraciones}...")
    
    # ETKDGv2 como método base
    df_etkdgv2 = pd.DataFrame(calculador_3d.pandas(generar_conformacion_etkdg(moleculas_para_conformaciones, 2)))
    resultados_descriptores[iteracion] = df_etkdgv2
    
    # Añadir resultados de ETKDGv1
    df_etkdgv1 = pd.DataFrame(calculador_3d.pandas(generar_conformacion_etkdg(moleculas_para_conformaciones, 1)))
    resultados_descriptores[iteracion] = pd.concat([resultados_descriptores[iteracion], 
                                                   df_etkdgv1.add_suffix("_ETKDG")], axis=1)
    
    # Añadir resultados de ETDG
    df_etdg = pd.DataFrame(calculador_3d.pandas(generar_conformacion_etdg(moleculas_para_conformaciones)))
    resultados_descriptores[iteracion] = pd.concat([resultados_descriptores[iteracion], 
                                                   df_etdg.add_suffix("_ETDG")], axis=1)
    
    # Añadir resultados de KDG
    df_kdg = pd.DataFrame(calculador_3d.pandas(generar_conformacion_kdg(moleculas_para_conformaciones)))
    resultados_descriptores[iteracion] = pd.concat([resultados_descriptores[iteracion], 
                                                   df_kdg.add_suffix("_KDG")], axis=1)
    
    # Añadir resultados de UFF
    df_uff = pd.DataFrame(calculador_3d.pandas(generar_conformacion_uff(moleculas_para_conformaciones)))
    resultados_descriptores[iteracion] = pd.concat([resultados_descriptores[iteracion], 
                                                   df_uff.add_suffix("_DGuff")], axis=1)
    
    # Añadir resultados de MMFF
    df_mmff = pd.DataFrame(calculador_3d.pandas(generar_conformacion_mmff(moleculas_para_conformaciones)))
    resultados_descriptores[iteracion] = pd.concat([resultados_descriptores[iteracion], 
                                                   df_mmff.add_suffix("_DGmmff")], axis=1)

print("Cálculo de descriptores completado.")

# Mostrar ejemplo de resultados
print(f"\nEjemplo de resultados (iteración 1): {resultados_descriptores[1].shape}")

# ----------------------------------------------------------------------------
# ANÁLISIS ESTADÍSTICO DE DISTRIBUCIONES
# ----------------------------------------------------------------------------

print("\n=== INICIANDO ANÁLISIS ESTADÍSTICO ===")

# Diccionario para almacenar estadísticas de distribución por molécula
estadisticas_distribucion = {}

# Calcular estadísticas para cada molécula
for indice_molecula in indices_moleculas:
    print(f"Analizando distribución para molécula {indice_molecula + 1}/{len(indices_moleculas)}...")
    
    # Recopilar valores de todos los descriptores para esta molécula en todas las iteraciones
    valores_molecula = pd.DataFrame()
    
    for iteracion in lista_iteraciones:
        serie_valores = pd.Series(resultados_descriptores[iteracion].iloc[indice_molecula])
        valores_molecula = valores_molecula.append(serie_valores, ignore_index=True)
    
    # Calcular estadísticas descriptivas (desde percentil 25 en adelante)
    estadisticas_distribucion[indice_molecula] = valores_molecula.describe()[3:]

print("Análisis estadístico completado.")

# Ejemplo de estadísticas para la molécula 2
print(f"\nEjemplo de estadísticas (molécula 2): {estadisticas_distribucion[2].shape}")

# ----------------------------------------------------------------------------
# PREPARACIÓN DE MATRIZ FINAL DE CARACTERÍSTICAS
# ----------------------------------------------------------------------------

print("\n=== PREPARANDO MATRIZ FINAL DE CARACTERÍSTICAS ===")

# Obtener estructura de nombres para las columnas finales
muestra_estadisticas = estadisticas_distribucion[1]
indices_estadisticas = list(range(len(muestra_estadisticas)))
indices_descriptores = list(range(len(muestra_estadisticas.T)))

# Generar nombres de columnas combinando estadística + descriptor
nombres_columnas_finales = []
for indice_descriptor in indices_descriptores:
    for indice_estadistica in indices_estadisticas:
        nombre_columna = (muestra_estadisticas.index[indice_estadistica] + 
                         muestra_estadisticas.columns[indice_descriptor])
        nombres_columnas_finales.append(nombre_columna)

print(f"Nombres de columnas generados: {len(nombres_columnas_finales)}")

# Construir matriz final de características estadísticas 3D
matriz_descriptores_estadisticos = pd.DataFrame()

for indice_molecula in indices_moleculas:
    # Recopilar todos los valores para esta molécula
    valores_molecula = pd.DataFrame()
    
    for iteracion in lista_iteraciones:
        serie_valores = pd.Series(resultados_descriptores[iteracion].iloc[indice_molecula])
        valores_molecula = valores_molecula.append(serie_valores, ignore_index=True)
    
    # Calcular estadísticas y reorganizar en formato vectorial
    estadisticas_molecula = valores_molecula.describe()[3:]
    fila_estadisticas = pd.DataFrame(estadisticas_molecula.values.reshape(1, 
                                   len(estadisticas_molecula) * len(estadisticas_molecula.T), 
                                   order='F'))
    
    matriz_descriptores_estadisticos = matriz_descriptores_estadisticos.append(fila_estadisticas, 
                                                                              ignore_index=True)

# Asignar nombres de columnas
matriz_descriptores_estadisticos.columns = [nombres_columnas_finales]

print(f"Matriz de descriptores estadísticos 3D: {matriz_descriptores_estadisticos.shape}")

# ----------------------------------------------------------------------------
# CÁLCULO DE DESCRIPTORES 2D COMPLEMENTARIOS
# ----------------------------------------------------------------------------

print("\n=== CALCULANDO DESCRIPTORES 2D COMPLEMENTARIOS ===")

# Usar moléculas originales sin conformación 3D
moleculas_2d = moleculas_para_conformaciones

# Crear calculador para descriptores 2D
calculador_2d = Calculator(descriptors, ignore_3D=True)

# Calcular descriptores 2D
matriz_descriptores_2d = calculador_2d.pandas(moleculas_2d)

print(f"Descriptores 2D calculados: {matriz_descriptores_2d.shape}")

# ----------------------------------------------------------------------------
# COMBINACIÓN Y LIMPIEZA DE DATOS FINALES
# ----------------------------------------------------------------------------

print("\n=== COMBINANDO Y LIMPIANDO DATOS FINALES ===")

# Combinar descriptores 2D y estadísticas 3D
matriz_completa_descriptores = pd.concat([matriz_descriptores_2d, matriz_descriptores_estadisticos], axis=1)

print(f"Matriz completa de descriptores: {matriz_completa_descriptores.shape}")

# Convertir a string para identificar valores no numéricos
matriz_limpieza = matriz_completa_descriptores.astype(str)

# Identificar celdas con texto (valores problemáticos)
mascara_texto = matriz_limpieza.apply(lambda columna: columna.str.contains('[a-zA-Z]', na=True))

# Eliminar filas/columnas con valores no numéricos
matriz_numerica = matriz_limpieza[~mascara_texto]

# Convertir de vuelta a numérico
matriz_numerica = matriz_numerica.astype(float)

# Rellenar valores faltantes con "NA"
matriz_numerica = matriz_numerica.fillna("NA")

print(f"Matriz después de limpieza: {matriz_numerica.shape}")

# ----------------------------------------------------------------------------
# COMBINACIÓN CON DATOS ORIGINALES Y EXPORTACIÓN
# ----------------------------------------------------------------------------

print("\n=== EXPORTANDO RESULTADOS FINALES ===")

# Combinar con el dataset original
dataset_final = pd.concat([df_original, matriz_numerica], axis=1)

# Exportar a archivo Excel
nombre_archivo_salida = 'analisis_descriptores_moleculares_completo.xlsx'
dataset_final.to_excel(nombre_archivo_salida, index=False)

print(f"Análisis completado exitosamente!")
print(f"Archivo exportado: {nombre_archivo_salida}")
print(f"Dataset final: {dataset_final.shape}")
print(f"- Moléculas procesadas: {len(indices_moleculas)}")
print(f"- Métodos de conformación utilizados: 6 (ETKDGv1, ETKDGv2, ETDG, KDG, UFF, MMFF)")
print(f"- Iteraciones por método: {numero_iteraciones}")
print(f"- Descriptores 3D base: {len(descriptores_filtrados)}")
print(f"- Total de características finales: {dataset_final.shape[1]}")

print("\n=== ANÁLISIS FINALIZADO ===")