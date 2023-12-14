#!/bin/bash

# Inserisci in valori_rispota la varietà di dimensionalità della matrice     
valori_risposta=(3 4 5 6 7 8 9 10 15 20 25 30 35 40 45 50 60 70 80 90 100)
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
    risposta_0="Problem1"
    risposta_1="no"
    risposta_2="$risposta_attuale"
    risposta_3="vandermode"
    echo -e "\033[1;33mCiclo $cont - Esecuzione con parametri : $risposta_2 e $risposta_3\033[0m"
    # Controlla se la compilazione è riuscita
    if [ $? -eq 0 ]; then
        echo -e "\033[;33mCompilazione completata con successo al ciclo $cont\033[0m"
        # Esegui l'eseguibile
expect <<EOF
        spawn ./main 
        expect "Prima domanda:"
        send "$risposta_0\r"
        expect "Seconda domanda:"
        send "$risposta_1\r"
        expect "Terza domanda:"
        send "$risposta_2\r"
        expect "Quarta domanda:"
        send "$risposta_3\r"
        expect eof 
EOF
    else 
        echo -e "\033[1;31mErrore durante la compilazione al ciclo $cont\033[0m"
    fi    

    # Appendere il risultato in base alla dimensionalità
    cp lufact.txt lufact_${risposta_3}_${risposta_attuale}.txt 
    mv /mnt/c/Users/Peppe/Desktop/AnalisiNumerica/lufact_${risposta_3}_${risposta_attuale}.txt /mnt/c/Users/Peppe/Desktop/AnalisiNumerica/$risposta_3/
    rm main
done

