#!/bin/bash

# Inizializza l'array di valori_risposta con numeri da 3 a 60
valori_risposta=()
for ((i=3; i<=60; i++)); do
    valori_risposta+=("$i")
done

# Contatore
cont=0
# Ciclo per ogni valore nella lista dei valori della risposta
for risposta_attuale in "${valori_risposta[@]}"; do
    # Incremento del contatore
    ((cont++))
    # Elenco dei file .f95 nella directory corrente
    file_list=$(ls *.f95)
    # Compila tutti i file .f95
    gfortran -fcheck=all -Wall -g -fbacktrace -o main $file_list -llapack -lblas
    risposta_1="Problem2"
    risposta_2="$risposta_attuale"
    echo -e "\033[1;33mCiclo $cont - Esecuzione con parametro : $risposta_2 \033[0m"
    # Controlla se la compilazione è riuscita
    if [ $? -eq 0 ]; then
        echo -e "\033[;33mCompilazione completata con successo al ciclo $cont\033[0m"
        # Esegui l'eseguibile
expect <<EOF
        spawn ./main 
        expect "Prima domanda:"
        send "$risposta_1\r"
        expect "Seconda domanda:"
        send "$risposta_2\r"
        expect eof 
EOF
    else 
        echo -e "\033[1;31mErrore durante la compilazione al ciclo $cont\033[0m"
    fi    

    # Appendere il risultato in base alla dimensionalità
    cp wilkin.txt wilkin_${risposta_attuale}.txt 
    cp lufact.txt lufact_wilkin_${risposta_attuale}.txt
    mv /mnt/c/Users/Peppe/Desktop/AnalisiNumerica/wilkin_${risposta_attuale}.txt /mnt/c/Users/Peppe/Desktop/AnalisiNumerica/wilkin/
    mv /mnt/c/Users/Peppe/Desktop/AnalisiNumerica/lufact_wilkin_${risposta_attuale}.txt /mnt/c/Users/Peppe/Desktop/AnalisiNumerica/wilkin/analisi_lufact/
    rm main
done

