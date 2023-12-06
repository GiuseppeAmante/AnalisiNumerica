#!/bin/bash

# Elenco dei file .f95 nella directory corrente
file_list=$(ls *.f95)

# Compila tutti i file .f95
#gfortran -o main $file_list
gfortran -fcheck=all -Wall -g -fbacktrace -o main $file_list
# Controlla se la compilazione è riuscita
if [ $? -eq 0 ]; then
  echo "Compilazione completata con successo."
  # Esegui l'eseguibile
  ./main
  vim factLU.txt

else
  echo "Errore durante la compilazione."
fi

