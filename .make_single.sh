#!/bin/bash

# Elenco dei file .f95 nella directory corrente
file_list=$(ls *.f95)
# Compila tutti i file .f95
gfortran -fcheck=all -Wall -g -fbacktrace -o main $file_list -llapack -lblas

if [ $? -eq 0 ]; then
    echo -e "\033[;33mCompilazione completata con successo al ciclo \033[0m"
    # Esegui l'eseguibile
    ./main 
else
    echo -e "\033[1;31mErrore durante la compilazione al ciclo \033[0m"
fi    
        
