Se añade el txt donde aparece el código usado para obtener la calidad por el programa de Trimmomatic.
Se procesan los datos a través del programa de Trimmomatic que lo que hará será eliminar adaptores que se han utilizado para poder secuenciar la muestra, recorte de bases de baja calidad, filtra lecturas cortas o de mala calidad y mantiene la coherencia entre pares, generando el archivo paired-end (utilizará para salmon) y separa las lecturas que han quedado impares en unpaired-end.
Para poder ver la calidad, realizamos FastQC en todas las muestras y MultiQC que nos realiza un resumen de la calidad de todas las muestras en un único documento htlm.
Solo se añade los resultados en MultiQC ya que se recogen todas las muestras y su calidad.
