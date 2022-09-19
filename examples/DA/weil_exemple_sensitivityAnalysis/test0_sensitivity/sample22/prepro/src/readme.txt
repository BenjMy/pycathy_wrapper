*******************************************************************************
*                                                                             *
*                             HAP versione 11.7                               *
*                                                                             *
*******************************************************************************

What's new in v11.7:

1. Nelle versioni precedenti si sono evidenziati alcuni problemi numerici dovuti 
all'arrotondamento dell'ultima cifra significativa in modo differente nei 
diversi ambienti quando la somma delle deviazioni è uguale nelle due direzioni; 
tali differenze però si ripercuotevano in modo macrospcopco nei risultati.
Si è scelto quindi di trattare il caso in cui le variabili a1 ed a2 sono uguali
in modo separato e definendo l'ugualgianza ponendo la differenza del valore 
assoluto minore di una certa soglia.




What's new in v11.5.1/v11.6:

1. Di alcuni prgrammi e' stato creato un main program che si chiama come nelle 
versioni precedenti che ha la sola funzione di chiamare una subroutine che 
contiene il listato del programma. 
In questo modo e' possbile compilare programmi in modo singolo come avveniva 
fino alla versione 11.4.1, oppure compilare un solo programma cppp.exe che 
permette di effettuare tutte le procedure necessarie per la definizione dei 
files di input di CATHY lanciando un solo eseguibile come avviene nella 
versione v11.5.
Cio' permette al bisogno di poter modificare la subroutine e poi poter avere le
modifche disponibili sia nel eseguibile unico sia in quelli singoli.

2. Introduzione nella subroutine wbb_sr.f90, di una procedura automatica per la
costruzione della gronda perimetrale del bacino come è stato implementato nella
versione 11.4.1.
E' possibile quindi definire automaticamente un percorso perimetrale di celle,
che può essere di forma regolare o irregolare, in modo tale da poter analizzare
bacini senza che siano ritagliati preventivamente.
La function catchment_cells e' stata traformata in subroutine; ciò ha permesso
di poter definire la cella di quota più bassa del bacino e in base ad essa in
modo automatico, in funzione di coefficienti definiti dall'utente nel file 
hap.in, di definire le quote da associare alle celle appartenenti alla gronda
e alla cella di chiusura.
Una piccola verifica poi permette di segnalare se i coefficienti sono stati
scelti opportunamente oppure se tali scelte potrebbero, nella fase di depit,che
interessa sicuramente le celle di gronda, delle modifiche ingiustificate anche
alla celle non appartenenti alla gronda.

3. La subroutine dsf.f90 è stata aggiornata in modo tale che possano essere 
determinati i numeri di Horton e di Strahler per ogni tratto di canale.
Sono state aggiunte tre variabili nel modulo mbbio e due record nella scritture
del file river_net.shp per immagazzinare i dati calcolati.
La definizione di queste due variabili per ciascuna cella viene effettuata 
solamente se il parametro ndcf (Non Dispersive Channel Flow) in hap.in è posto 
uguale ad 1. Nonostante i metodi dispersivi siano selezionati operando su altri 
parametri, è bene mantenere sempre l'impostazione di ndcf pari ad 1. Si consiglia
di cambiarla secondo le esigenze solamente se si utilizzano metodi dispersivi o 
parzialmente dispersivi.

4. Il programma rbb ora offre la possibilità, quando si sta utilizzando un
metodo non dispersivo, di poter scrivere un file unico per quelle variabili 
che per ogni cella sono generalmente immagazzinate secondo le direzioni 
cardinali e quelle diagonali.

5. Il file batch compf90.bat permette da schermo di scegliere quale programma 
si vuole compilare e lancia la stringa per la compilazione con il compilatore
lf95.

*******************************************************************************

