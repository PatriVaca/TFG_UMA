#######################################################################
##      Escrito por Patricia del Carmen Vaca Rubio, 2023             ##
##  Estudiante del grado de Bioquímica con mención en Biotecnología  ##
##          Universidad de Málaga, Andalucía, España                 ##
#######################################################################

####Sección 1####

BiocManager::install(version = "3.16")
BiocManager::install("biomaRt")
BiocManager::install("tibble")

library(tibble)
require(biomaRt)

#Primero, se importaron en R los datos de los archivos .txt (matrices de conteo) con 'readLines'
matriz_conteo_Gg=readLines("matrix_count_Gg.txt") #Matriz de conteo del pollo (Gallus gallus)
matriz_conteo_Mm=readLines("matrix_count_Mm.txt") #Matriz de conteo del ratón (Mus musculus)

#Segundo, como el output de readLines en cada caso fue un vector de caracteres de una sola columna, se obtuvo en su lugar cada matriz de conteo como un data.frame con read.delim(), puesto que esta función separa el vector por tabulación (\t) en columnas
genes_Gg <- read.delim("matrix_count_Gg.txt", header=TRUE) #header=TRUE para mantener los nombres originales de las columnas
genes_Mm <- read.delim("matrix_count_Mm.txt", header=TRUE)

names(genes_Gg)[1]="ENSEMBL_ID" #para cambiar el nombre de la primera columna de Geneid a ENSEMBL_ID
names(genes_Mm)[1]="ENSEMBL_ID"

#Luego, se realizó una conexión a la base de datos de BioMart, Ensembl, mediante useMart --> para recuperar toda la información (con dataset) de esa base de datos (Ensembl) de los genes de Mus musculus.
ensembl_Mm <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

#A continuación, se recuperaron en el dataframe annot_Mm los atributos "ensembl_gene_id" y "external_gene_name" de la base de datos Ensemble (a partir de la búsqueda de antes, mart=ensembl_Mm) equivalentes a los genes de la columna ENSEMBL_ID del dataframe de genes_Mm
annot_Mm <- getBM(
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name'),
  filters = 'ensembl_gene_id',
  values = genes_Mm$ENSEMBL_ID,
  mart = ensembl_Mm)

#Posteriormente, en el caso del pollo, la versión de la anotación original (resultado del RNA-seq) de la matriz de conteo se correspondía con la GRCg6a (que se encuentra en la versión de Ensembl 106) 
#Sin embargo, a la versión de Ensembl a la que se accede por defecto a través de useMart es a la más reciente (en este caso,fue a la versión 109, que tenía la anotación del genoma de pollo GRCg7b)
#Por lo tanto, para la correcta anotación en R del dataframe de pollo, se buscó por internet a qué versión de Ensembl o "Ensemble release version" pertenecía mi anotación deseada, la GRCg6a.
#Con el comando "listEnsemblArchives()" se pudieron ver el nombre, la fecha, el URL y la versión de cada "Ensembl release version". Se copió el URL de la versión de interés (en este caso, la Ensembl release 106)
listEnsemblArchives()
listMarts() #Para ver a qué bases de datos de BioMart se podía acceder mediante el paquete "biomaRt".

#La base de datos de Ensembl (de las 4 que se mostraban como output de listMarts()) a la que accedí a través de biomaRt fue a la de "ENSEMBL_MART_ENSEMBL", ya que era la que contenía la información de los genes
#A continuación, se creó un objeto mart pero, esta vez, sin usar useMart (ya que este comando no presentaba la opción de recuperar versiones anteriores de Ensembl). En su lugar, se utilizó useEnsembl:

ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://apr2022.archive.ensembl.org") #La clave para recuperar la Ensembl release 106 consistió en poner host=URL que copiamos del output de listEnsemblArchives().
ensembl_Gg <- useDataset(dataset = "ggallus_gene_ensembl", mart = ensembl) #Para seleccionar el dataset deseado (de pollo) de la base de datos de BioMart seleccionada con el comando anterior.

annot_Gg <- getBM(
  attributes = c(
    'ensembl_gene_id',
    'external_gene_name'),
  filters = 'ensembl_gene_id',
  values = genes_Gg$ENSEMBL_ID,
  mart = ensembl_Gg) #Se recuperó la anotación deseada para el pollo

#Seguidamente, se unieron los dos data.frames, genes_Mm y annot_Mm, en uno por los valores de la columna "ensembl_gene_id" (desapareciendo este encabezado) y manteniendo todas las filas (all.x=T) de la columna ENSEMBL_ID.
merged_Mm <- merge(
  x = genes_Mm,
  y =  annot_Mm,
  by.y = 'ensembl_gene_id',
  all.x = T,
  by.x = 'ENSEMBL_ID') 

#Se cargó el paquete dplyr para usar la función relocate()
library(dplyr) # from version 1.0.0 

merged_Mm2 <- merged_Mm %>%
  relocate(external_gene_name) #Se colocó la columna de external_gene_name a la izquierda del todo (como primera columna)

merged_Mm2 <- select (merged_Mm, -ENSEMBL_ID) #Se quitó la columna ENSEMBL_ID

#A continuación, se hizo lo mismo con el pollo:

merged_Gg <- merge(
  x = genes_Gg,
  y =  annot_Gg,
  by.y = 'ensembl_gene_id',
  all.x = T,
  by.x = 'ENSEMBL_ID') 

merged_Gg2 <- merged_Gg %>%
  relocate(external_gene_name)
merged_Gg2 <- select (merged_Gg2, -ENSEMBL_ID)

#Se quitaron los caracteres en blanco, es decir, los genes no anotados de cada especie:
merged_Gg3 <- merged_Gg2[merged_Gg2$external_gene_name != '',]
merged_Mm3 <- merged_Mm2[merged_Mm2$external_gene_name != '',]

#Para ver qué genes eran ortólogos entre ratón y pollo:
common_external_gene_name <- intersect(merged_Mm3$external_gene_name, merged_Gg3$external_gene_name)

#Se unieron en un solo dataframe únicamente las filas del dataframe de pollo y de ratón que fueran compartidas (compartían el nombre recuperado con external_gene_name, es decir, eran ortólogos) con inner_join()
merged_orthologues <- inner_join(merged_Mm2, merged_Gg2, by = "external_gene_name") #Con genes no anotados
merged_annotated_orthologues <- inner_join(merged_Mm3, merged_Gg3, by = "external_gene_name") #Sin genes no anotados

#Sin embargo, para encontrar los ortólogos correctamente, se necesitaba que la anotación de los gene_external_name fuera igual en pollo y en ratón. Para ello, se pasaron todos los external_gene_name a minúsculas con tolower()
merged_Gg3$external_gene_name=tolower(merged_Gg3$external_gene_name)
merged_Mm3$external_gene_name=tolower(merged_Mm3$external_gene_name)

merged_annotated_orthologues_2 <- select(merged_annotated_orthologues, c(external_gene_name,Mm_PE_1:Mm_PE_3,PE_1:PE_3)) #Para quedarnos únicamente con las columnas de las tres réplicas biológicas del proepicardio (PE) de ratón y pollo en la matriz final

write.csv(merged_annotated_orthologues_2, "C:/Users/patri/Desktop/UMA 2/4BQ/TFG/conversion_TFG/matriz_conteo_TFG//merged_annotated_orthologues_2.csv") #Se guardó el dataframe como un archivo .csv

####Sección 2####

#A continuación, se cambiaron los external_gene_names por el ENSEMBLID de humano correspondiente, para poder utilizar luego más herramientas bioinformáticas (basadas en humano)

ensembl_Hs <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl') #Uso del dataset "hsapiens_gene_ensembl" que se correspondía con la versión de anotación "GRCh38.p13"

annot_Mm_Hs <- getBM(
  attributes = c('ensembl_gene_id',
                 'external_gene_name',
                 'hsapiens_homolog_ensembl_gene'),
  filters = 'ensembl_gene_id',
  values = genes_Mm$ENSEMBL_ID,
  mart = ensembl_Mm) #Para recuperar la anotación de los ortólogos en humano para los genes de ratón (también podría ser de pollo)

#Posteriormente, se unió el dataframe de la anotación de humano con el dataframe obtenido anteriormente "merged_annotated_orthologues_2"

annot_Mm_Hs$external_gene_name<-tolower(annot_Mm_Hs$external_gene_name)
annot_Hs<-select(annot_Mm_Hs, -ensembl_gene_id)

annot_Hs <- annot_Hs[!duplicated(annot_Hs$external_gene_name), ]

merged_annotated_orthologues_3 <- inner_join(annot_Hs, merged_annotated_orthologues_2, by = "external_gene_name")
merged_annotated_orthologues_3 <- select (merged_annotated_orthologues_3, -external_gene_name)

#Sin embargo, se observó que merged_annotated_orthologues_3 tenía genes sin anotar al no recuperarse una correspondencia entre algunos "external gene names" y su ENSEMBL ID humano. 

merged_annotated_orthologues_4 <- merged_annotated_orthologues_3[merged_annotated_orthologues_3$hsapiens_homolog_ensembl_gene != '',] #Para eliminar los genes sin anotación de la matriz final de conteo

#Se guardó la matriz de conteo finalmente procesada en formato .csv para poder cargarlo posteriormente en la herramienta web iDEP.96
write.csv(merged_annotated_orthologues_4, "C:/Users/patri/Desktop/UMA 2/4BQ/TFG/conversion_TFG/matriz_conteo_TFG//merged_annotated_orthologues_Hs_anotados.csv")


####Sección 3####

#Tratamiento de los datos obtenidos con iDEP.96 con un FDRcutoff=0,05 para poder cargarlos posteriormente en otra herramienta web, ShinyGO

#Se obtuvo, por un lado, un archivo que solo que solo tenga los genes downregulated que aparecen en iDEP.96, a partir del archivo All_gene_lists_GMT.

setwd("C:/Users/patri/Desktop/UMA 2/4BQ/TFG/iDEP.96/Con ENSEMBLID humano/FDR cutoff 0,05") #Para cambiar al nuevo directorio de trabajo

matriz_shinyGO_0.05=readLines("All_gene_lists_GMT.txt") #Se importó en R el archivo de texto generado por iDEP.96 que contenía todos los genes expresados diferencialmente (entre el PE de pollo y de ratón) 
matriz_shinyGO_0.05=read.delim("All_gene_lists_GMT.txt", header=FALSE) #header=FALSE para no mantener los nombres originales de las columnas

matriz_shinyGO_upregulated_0.05=matriz_shinyGO_0.05[-2,] #Se obtuvo una matriz que solo contenía los genes regulados al alza
matriz_shinyGO_downregulated_0.05=matriz_shinyGO_0.05[-1,] #Se obtuvo una matriz que solo contenía los genes regulado a la baja

#Se exportaron las matrices obtenidas en formato .csv
write.csv(matriz_shinyGO_upregulated_0.05, "C:/Users/patri/Desktop/UMA 2/4BQ/TFG/iDEP.96/Con ENSEMBLID humano/FDR cutoff 0,05//todoslosIDs_upregulated_0,05.csv")
write.csv(matriz_shinyGO_downregulated_0.05, "C:/Users/patri/Desktop/UMA 2/4BQ/TFG/iDEP.96/Con ENSEMBLID humano/FDR cutoff 0,05//todoslosIDs_downregulated_0,05.csv")
write.csv(matriz_shinyGO_0.05, "C:/Users/patri/Desktop/UMA 2/4BQ/TFG/iDEP.96/Con ENSEMBLID humano/FDR cutoff 0,05//todoslosIDs_DEG_0,05.csv")




