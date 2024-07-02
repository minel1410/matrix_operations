Matrix operations
Minel Salihagić

Tema projekta je implementacija klase Matrica koja čuva matricu proizvoljnih dimenzija. 
Naglasak je bio na brzom množenju matrice, determinanti i operatoru '>>'.
S obzirom da nije naglašeno da li klasa treba biti generička, ja sam se odlučio za implementaciju 
matrice preko vektora vektora, koji čuvaju doubleove. Pošto je u pitanju double podržana je i 
automatska konverzija iz inta u double. Ono što je bitno naglasiti je da klasa sadrži član 'skalar' 
koji čuva informaciju o tome da li je matrica stvarno matrica ili je broj predstavljen u obliku 
matrice, što je bilo potrebno za implementaciju operatora '>>', koji je i evaluator izraza sa 
matricama.

Postoji i par pomoćnih privatnih funkcija koji služe za lakše implementiranje ostalih funkcija. 
Većina imena funkcija kao i njihova implementacija su jednostavne pa ih neću detaljnije
pojašnjavati ovdje, kratko objašnjenje je dato u komentarima za svaku funkciju. 
Funkcija 'Strassen' koja prima dvije matrice i vraća matricu je funkcija koja vrši množenje dvije 
matrice u vremenu 𝑂(𝑛
2,81) što je malo brže u odnosu na standardni algoritam za množenje 
𝑂(𝑛
3
). 
Funkcija je implementirana metodom „podijeli pa vladaj“ kao što je traženo u zadatku. 
Radi tako što podjeli matrice na po 4 bloka, te svaki blok množi zasebno. Ako je blok veći od 2x2, 
što je i bazni slučaj, matrica se rekurzivno poziva na te blokove. Nakon završetka algoritma, 
rezultat su 4 matrice, koje su podmatrice rezultujuće. Potrebno je samo spojiti ih u jednu. S 
obzirom da Strassenov algoritam radi samo za kvadratne matrice, čija je dužina stepen broja 2, 
potrebno je dodati određen broj 0-kolona i redova da bi dobili takve matrice. Za to služe funkcije 
'dodajRedNaKraj ()' i ostale. Funkcija 'obrisiVisak ()' miče redove i kolone koji su višak nakon 
množenja. Funkcija 'ispraviNule ()' služi za zaokruživanje veoma malih brojeva na 0 koji nastaju 
zbog floating point aritmetike prilikom množenja. 
Sve funkcije su implementirane na način da rade i sa konstantnim i nekonstantim matricama, 
odnosno da rade bez promjene podataka.

Što se tiče konstruktora, postoji konstruktor bez parametara, konstruktor sa 2 parametra koji 
pravi jediničnu matricu proslijeđenih dimenzija te baca izuzetak ukoliko je parametar negativan 
ili nula, konstruktor kopije koji prima skalar, te ga čuva kao matricu. Zatim konstruktor koji prima 
matricu te dva parametra, i vraća istu tu matricu bez reda i kolone koji su proslijeđeni. Služi kod 
računanja determinante. Pored toga implementirano je „The Big Five“, te također operator 
dodjele za skalarne matrice.

Funkcija 'Inverzna ()' vraća matricu koja je inverz proslijeđene matrice. Prvo provjeri da li je 
matrica invertibilna računajući determinantu, te ako nije baca izuzetak. Radi po formuli:
𝐴
−1 =
1
𝑑𝑒𝑡𝐴 ∙ 𝑎𝑑𝑗𝐴
gdje je 𝑎𝑑𝑗𝐴 transponovana matrica algebarskih komplemenata, odnosno matrica čiji je svaki 
element jednak vrijednosti determinante matrice bez kolone i reda u kojoj se nalazi taj element.
Funkcija radi u 𝑂(𝑛
3
) za n x n matricu.
Funkcija za računanje determinante je implementirana da radi Gausovom metodom sa 
pivotiranjem, gdje se matrica svodi na gornju trougaonu te radi u vremenu 𝑂(𝑛
3
).

S obzirom da je potrebna provjera da li je skalar true ili false u nekim funkcijama koje nisu 
deklarisane kao prijatelji klase, napravljen je geter za skalar. 
Implementirani su i operatori „ ==, !=, +=, -=, *=, /=, +, -, *, /, ^“ koji rade i za slučajeve kada je 
jedna od proslijeđenih matrica skalar. Operator '^' radi po algoritmu brzog stepenovanja. Ako je 
eskponent manji od nula, matrica se prvo invertuje, a eksponent se množi sa -1 te se poziva 
operator, sada na invertovanu matricu i pozitivan eksponent. Ako je eksponent 0, vraća se 
jedinična matrica, a ako je veći od nula posmatra se binarni zapis eksponenta, odnosno rezultat 
pri dijeljenju sa 2. Ako je broj djeljiv sa 2, odnosno ako je u binarnom zapisu cifra 2, matrica se 
množi baznom matricom, odnosno onom matricom koja je proslijeđena te se eksponent podijeli 
sa 2, i tako sve dok je eksponent veći od 0. 

Osim toga podržani su i operatori za ispis i upis matrica. Zadnja funkcija je 'uzmiMatricu ()' koja 
je pomoćna funkcija kod operatora '>>', koja služi za izdvajanje samo matrice iz ulaznog toka, 
kada naiđe na '['. Ovdje je bitno napomenuti da se kod unošenja izraza podrazumijeva da znak '[' 
označava početak unosa matrice, te će javiti grešku ako matrica nije unesena u ispravnom 
formatu.

Operator '>>' radi po „shunting yard“ algoritmu, na sličan način na koji je implementiran i 
evaluator na vježbama, te su implementirane pomoćne funkcije koje olakšavaju implementaciju 
operatora. 
U main funkciji su napravljene funkcije koje služe za generisanje random intova i doubleova, te
jedna funkcija koja generiše random matricu. Napravljeno je par for petlji koje služe za testiranje 
nekih funkcija na matricama koje su random generisane. Također u komentaru je napisano i 5
izraza koji se mogu upisati da bi se testirala funkcionalnost operatora '>>'. Potrebno je 
napomenuti da kod operatora '>>', kod unošenja inverzne matrice potrebno je napisati u obliku 
„[1 2; 3 4] ^(-1)“.
Minel Salihagić, 5930/M
Strukture podataka i algoritmi, 17.1.2023.
