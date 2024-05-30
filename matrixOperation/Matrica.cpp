#ifndef MATRICA_CPP
#define MATRICA_CPP
#include "Matrica.h"
using namespace std;

// Funkcija koja provjerava da li je broj stepen od 2
bool stepenOdDva(int x) { return (x & (x - 1) && x > 0) == 0; }

bool Matrica::jeLiKvadratna() const {

    if (skalar)
        throw "Matrica je skalar.";

    return (broj_redova == broj_kolona);
}

// Matrica koja prosiruje proslijedjenu matricu do kvadratne matrice dodavajuci nula redove i kolone
Matrica Matrica::napraviKvadratnu() const {

    if (jeLiKvadratna())
        throw "Matrica je vec kvadratna";

    Matrica kvadratna = *this;
    if (broj_redova < broj_kolona) {

        vector<double>nule(broj_kolona, 0);
        while (kvadratna.broj_redova < kvadratna.broj_kolona)
            kvadratna.dodajRedNaKraj(nule);
    }
    else {
        vector<double>nule(broj_redova, 0);
        while (kvadratna.broj_kolona < kvadratna.broj_redova)
            kvadratna.dodajKolonuNaKraj(nule);
    }

    return kvadratna;

}

// Matrica koja se koristi kod Strassena da bi se proslijedjene matrice prosirile
Matrica Matrica::napraviStepenOdDva() const {

    Matrica ret = *this;
    if (!ret.jeLiKvadratna())
        ret = napraviKvadratnu();

    // Dodajemo nula red sve dok broj redova ne postane stepen od 2
    vector<double>nulaRed(ret.broj_kolona, 0);

    while (!stepenOdDva(ret.broj_redova))
        ret.dodajRedNaKraj(nulaRed);

    // Isto radimo i za kolone
    vector<double>nulaKolona(ret.broj_redova, 0);

    while (!stepenOdDva(ret.broj_kolona))
        ret.dodajKolonuNaKraj(nulaKolona);

    return ret;

}

// Funkcija koja se koristi za mnozenje kod Strassena, i uklanja redove ili kolone koje su nastale prosirenjem matrica
void Matrica::obrisiVisak(const Matrica& A, const Matrica& B) {

    while (broj_redova > A.broj_redova) {

        elementi.pop_back();
        broj_redova--;
    }

    while (broj_kolona > B.broj_kolona) {

        for (size_t i = 0; i < broj_redova; i++)
            elementi[i].pop_back();
        broj_kolona--;
    }
}


// Privatna funkcija koja mnozi dvije matrice koristeci Strassenov algoritam
Matrica Strassen(const Matrica& a, const Matrica& b) {

    // Ako matrice nisu odgovarajuce za mnozenje baca se izuzetak
    if (a.broj_kolona != b.broj_redova)
        throw "Matrice nisu odgovarajuce za mnozenje.";

    // Pravimo kvadratne matrice od matrica a i b, koje imaju broj redova koji je stepen od dva
    Matrica A = a.napraviStepenOdDva();
    Matrica B = b.napraviStepenOdDva();

    // Ako u nekom slucaju matrice nemaju isti broj redova nakon prosirenja, provjeravamo koju treba dodatno prosiriti
    if (A.broj_redova != B.broj_redova) {
        if(A.broj_redova < B.broj_redova) {

            while(A.broj_redova < B.broj_redova)
                A.dodajRedNaKraj(vector<double>(A.broj_kolona, 0));

            A = A.napraviStepenOdDva();
        }
        else {

        while(A.broj_redova > B.broj_redova)
                B.dodajRedNaKraj(vector<double>(B.broj_kolona, 0));

            B = B.napraviStepenOdDva();
        }
    }

    int K = A.broj_redova;

    // Pravimo povratnu matricu
    Matrica C(K, K);


    // Bazni slucajevi za koje se proizvod racuna u konstantnom vremenu
    if (K == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return C;
    }
    if (K == 2) {
        C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0];
        C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1];
        C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0];
        C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1];
    }
    int k = K/2;

    // Dijelimo matrice A i B na 4 jednake matrice
    Matrica A11(k, k);
    Matrica A12(k, k);
    Matrica A21(k, k);
    Matrica A22(k, k);

    Matrica B11(k, k);
    Matrica B12(k, k);
    Matrica B21(k, k);
    Matrica B22(k, k);


    // Dodajemo elemente u podmatrice
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {

            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j+k];
            A21[i][j] = A[i+k][j];
            A22[i][j] = A[i+k][j+k];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j+k];
            B21[i][j] = B[i+k][j];
            B22[i][j] = B[i+k][j+k];
        }
    }

        Matrica S1 = A12 - A22;
        Matrica S2 = B21 + B22;
        Matrica S3 = A11 + A22;
        Matrica S4 = B11 + B22;
        Matrica S5 = A11 - A21;
        Matrica S6 = B11 + B12;
        Matrica S7 = A11 + A12;
        Matrica S8 = B12 - B22;
        Matrica S9 = B21 - B11;
        Matrica S10 = A21 + A22;

        Matrica M1 = Strassen(S1, S2);
        Matrica M2 = Strassen(S3, S4);
        Matrica M3 = Strassen(S5, S6);
        Matrica M4 = Strassen(S7, B22);
        Matrica M5 = Strassen(A11, S8);
        Matrica M6 = Strassen(A22, S9);
        Matrica M7 = Strassen(S10, B11);


        // Racunamo podmatrice rezultata
        Matrica C11 = M1 + M2 - M4 + M6;
        Matrica C12 = M4 + M5;
        Matrica C21 = M6 + M7;
        Matrica C22 = M2 - M3 + M5 - M7;

        // Objedinjujemo podmatrice u jednu
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                C[i][j] = C11[i][j];
                C[i][j+k] = C12[i][j];
                C[i+k][j] = C21[i][j];
                C[i+k][j+k] = C22[i][j];
            }
        }
        return C;
}

// Privatna funkcija koja sluzi da veoma male brojeve koji zaokruzi na 0, zbog floating point aritmetike
void Matrica::ispraviNule() {

    for (size_t i = 0; i < broj_redova; i++)
        for (size_t j = 0; j < broj_kolona; j++)
            if (abs(0 - elementi[i][j]) < 0.000001)
                elementi[i][j] = 0;

}

// Konstruktor bez parametara koji pravi praznu matricu, sluzi za unosenje matrice preko operator  >>
Matrica::Matrica(): broj_redova(0), broj_kolona(0) {};

// Konstruktor sa dva parametra koji pravi jedinicnu matricu brRedova x brKolona
Matrica::Matrica(int brRedova, int brKolona) {

    if(brRedova <= 0 || brKolona <= 0)
        throw "Broj redova i kolona matrice mora biti pozitivan broj.";

    broj_kolona = brKolona;
    broj_redova = brRedova;

    for (size_t i = 0; i < broj_redova; i++) {
        vector<double>red;
        for (size_t j = 0; j < broj_kolona; j++) {
            if (i == j)
                red.push_back(1);
            else
                red.push_back(0);
        }
        elementi.push_back(red);
    }
}

// Konstruktor koji pretvara skalar u matricu 1x1
Matrica::Matrica(const double &skalar) : broj_redova(1), broj_kolona(1), skalar(true) {
    elementi.push_back(vector<double>({skalar}));
}

// Operator dodjele koji pretvara skalar u matricu 1x1
Matrica& Matrica::operator = (const double& skalar) {

    broj_redova = 1;
    broj_kolona = 1;
    this->skalar = true;
    elementi.clear();
    elementi.push_back(vector<double>({skalar}));
    return *this;

}

// Verzija konstruktora kopije koji vraca matricu bez reda i kolone sa indeksima I i J
// Koristi se samo kod racunanja determinante
Matrica::Matrica(const Matrica& m, int I, int J) {

    if(size_t(I) > m.broj_redova || size_t(J) > m.broj_kolona || I < 0 || J < 0)
        throw "Dati indeksi se ne nalaze u matricu";

    skalar = m.skalar;
    *this = m;
    for (size_t i = 0; i < broj_redova; i++)
        elementi[i].erase(elementi[i].begin() + J);

    elementi.erase(elementi.begin() + I);
    broj_redova--;
    broj_kolona--;

}

//  Konstruktor kopije
Matrica::Matrica(const Matrica& kopija) {

    if(this != &kopija) {

        broj_redova = kopija.broj_redova;
        broj_kolona = kopija.broj_kolona;
        skalar = kopija.skalar;
        elementi = kopija.elementi;

    }
}


// Pomjerajuci konstruktor kopije
Matrica::Matrica(Matrica&& kopija) {

    if(this != & kopija) {

        skalar = kopija.skalar;
        elementi = move(kopija.elementi);
        broj_redova = kopija.broj_redova;
        broj_kolona = kopija.broj_kolona;

        kopija.broj_redova = 0;
        kopija.broj_kolona = 0;
        kopija.skalar = false;
    }
}
// Operator dodjele
Matrica& Matrica::operator = (const Matrica& rhs) {

    if(this != &rhs) {

        skalar = rhs.skalar;
        broj_redova = rhs.broj_redova;
        broj_kolona = rhs.broj_kolona;
        elementi = rhs.elementi;

    }

    return *this;
}

// Pomjerajuci operator dodjele
Matrica& Matrica::operator = (Matrica&& rhs) {

    if (this != &rhs) {

        skalar = rhs.skalar;
        broj_redova = rhs.broj_redova;
        broj_kolona = rhs.broj_kolona;
        elementi = move(rhs.elementi);

        rhs.broj_redova = 0;
        rhs.broj_kolona = 0;
        rhs.skalar = false;
    }

    return *this;
}

// Destruktor za klasu
Matrica::~Matrica() {

    broj_redova = 0;
    broj_kolona = 0;

    // Brisemo sve elemente iz vektora
    elementi.clear();

    // S obzirom da je vektor prazan, smanjujemo mu velicinu
    elementi.shrink_to_fit();
}

// Funkcija koja vraca transponovanu matricu, u odnosu na proslijedjenu matricu
Matrica Matrica::Transponovana() const {

    if (skalar)
        throw "Matrica je skalar";

    Matrica m(broj_kolona, broj_redova);

    // Idemo po kolonama te dodajemo elemente u matricu
    for (size_t i = 0; i < broj_kolona; i++)
        for (size_t j = 0; j < broj_redova; j++)
            m[i][j] = elementi[j][i];
    return m;
}

// Funkcija za trazenje inverzne matrice
Matrica Matrica::Inverzna() const {

    if (skalar)
        throw "Matrica je skalar";

    // Funkcija radi po formuli A^(-1) = 1/detA * |A|^T, pa je potrebno naci determinantu
    double det = Determinanta();

    // Funkcija je invertibilna ako je determinanta razlicita od 0
    if (det == 0)
        throw "Matrica nije invertibilna";

    det = 1.0 / det;
    Matrica inverzna(broj_redova, broj_kolona);

    // Racunamo algebarski komplement za svaki element matrice
    for (size_t i = 0; i < broj_redova; i++) {
        for (size_t j = 0; j < broj_kolona; j++) {

            double trenutnaDet = Matrica(*this, i, j).Determinanta();

            // Mnozimo trenutnu determinantu sa 1 ili -1 u ovisnosti od parnosti indeksa
            trenutnaDet *= ((i+j) % 2 == 0 ? 1 : -1);
            inverzna[i][j] = trenutnaDet;
        }
    }

    // Na kraju transponujemo matricu, te mnozim sa det = 1 / det
    Matrica inverznaT = inverzna.Transponovana();
    inverznaT = det * inverznaT;

    for (size_t i = 0; i < inverznaT.broj_redova; i++) {
        for (size_t j = 0; j < inverznaT.broj_kolona; j++) {
            if (inverznaT[i][j] == -0)
                inverznaT[i][j] = 0;
        }

    }
    return inverznaT;
}

// Geter za skalar
bool Matrica::getSkalar() const { return skalar; }

// Funkcija koja dodaje proslijedjeni vektor kao red na kraj matrice
void Matrica::dodajRedNaKraj(vector<double>red) {

    if (skalar)
        throw "Matrica je skalar";

    // U slucaju da dati vektor nema duzinu kao i matrica, ne moze se dodati
    if (red.size() != broj_kolona)
        throw "Dati red se ne moze dodati u matricu.";

    elementi.push_back(red);
    broj_redova++;

}

// Funkcija koja dodaje proslijedjeni vektor kao red na pocetak matrice
void Matrica::dodajRedNaPocetak(vector<double>red) {

    if (skalar)
        throw "Matrica je skalar";

    // U slucaju da dati vektor nema duzinu kao i matrica, ne moze se dodati
    if (red.size() != broj_kolona)
        throw "Dati red se ne moze dodati u matricu.";

    elementi.insert(elementi.begin(), red);
    broj_redova++;
}

// Funkcija koja dodaje proslijedjeni vektor kao kolonu na kraj matrice
void Matrica::dodajKolonuNaKraj(vector<double>kolona) {

    if (skalar)
        throw "Matrica je skalar";

    // U slucaju da dati vektor nema duzinu kao i matrica, ne moze se dodati
    if (kolona.size() != broj_redova)
        throw "Data kolona se ne moze dodati u matricu.";

    for (size_t i = 0; i < broj_redova; i++)
        elementi[i].push_back(kolona[i]);
    broj_kolona++;
}

//  Funkcija koja dodaje proslijedjeni vektor kao kolonu na pocetak matrice
void Matrica::dodajKolonuNaPocetak(vector<double>kolona) {

    if (skalar)
        throw "Matrica je skalar";

    if (kolona.size() != broj_redova)
        throw "Data kolona se ne moze dodati u matricu.";

    for (size_t i = 0; i < broj_redova; i++)
        elementi[i].insert(elementi[i].begin(), kolona[i]);
    broj_kolona++;
}

//  Funkcija koja racuna determinantu matrice pomocu Gausove eliminaicije
double Matrica::Determinanta() const {

    if(skalar)
        throw "Matrica je skalar";

    //  Da li je matrica kvadratna
    if (!jeLiKvadratna()) {
        throw "Matrica nije kvadratna pa nema determinantu";
    }

    Matrica m = *this;

    double det = 1;
    for (size_t k = 0; k < broj_kolona; k++) {
        int pivotRed = k;
        //  Trazimo pivot red, tj. red sa najvecim elementom u trenutnoj koloni
        for (size_t i = k + 1; i < broj_redova; i++) {
            if (fabs(m[i][k]) > fabs(m[pivotRed][k])) {
                pivotRed = i;
            }
        }
        if (fabs(m[pivotRed][k]) < 1e-10) {
            return 0;
        }
        if (size_t(pivotRed) != k) {
            //  Mijenjamo trenutni red sa pivot redom
            for (size_t j = k; j < broj_kolona; j++) {
                swap(m[k][j], m[pivotRed][j]);
            }
            det *= -1;
        }
        det *= m[k][k];

        //   Oduzimamo pivot red pomonozen faktorom od svakod reda ispod
        for (size_t i = k + 1; i < broj_redova; i++) {
            double faktor = m[i][k] / m[k][k];
            for (size_t j = k + 1; j < broj_kolona; j++) {
                m[i][j] -= faktor * m[k][j];
            }
        }
    }
    return det;
}

// Operator koji vraca vektor iz matrice po referenci
vector<double>& Matrica::operator[] (int indeks) {

    if (size_t(indeks) >= broj_redova)
        throw "Indeks je van opsega matrice.";

    else if (indeks < 0)
        throw "Indeks ne moze biti negativan";

    return elementi[indeks];
}

// Operator koji vraca konstantan vektor iz konstante matrice
const vector<double>& Matrica::operator[] (int indeks) const {

    if (size_t(indeks) >= elementi.size())
        throw "Indeks je van opsega matrice.";

    else if (indeks < 0)
        throw "Indeks ne moze biti negativan";

    return elementi[indeks];
}

bool Matrica::operator == (const Matrica& m) {

    if (broj_redova != m.broj_redova || broj_kolona != m.broj_kolona)
        return 0;

    for (size_t i = 0; i < broj_redova; i++) {
        for (size_t j = 0; j < broj_kolona; j++)
            if (elementi[i][j] != m[i][j])
                return false;
    }
    return 1;
}

bool Matrica::operator != (const Matrica& m) { return !(*this == m); }

// Sabiranje date matrice sa drugom matricom
void Matrica::operator += (const Matrica& m) { *this = *this + m; }

// Oduzimanje matrice od date matrice
void Matrica::operator -= (const Matrica& m) { *this = *this - m; }

// Mnozenje matrice drugom matricom
void Matrica::operator *= (const Matrica& m) {

    // Ovdje nije potrebno provjeravati podudarnost jer ce se to izvrsiti u operator *
    *this = *this * m;
}


// Mnozenje matrice skalarom
void Matrica::operator *= (double skalar) { *this = *this * skalar; }

void Matrica::operator /= (double skalar) { *this = *this / skalar; }

// Operator + za sabiranje matrica
Matrica operator + (const Matrica& m1, const Matrica& m2) {

    // Provjeravamo da li su matrice A i B zapravo skalari u matricnoj formi
    if (m1.skalar && m2.skalar)
        return Matrica(m1[0][0] + m2[0][0]);

    // Provjera podudarnosti matrica
    if (m1.broj_redova != m2.broj_redova || m1.broj_kolona != m2.broj_kolona)
        throw "Nepodudarne dimenzije matrica za sabiranje!";

    Matrica rez(m1.broj_redova, m1.broj_kolona);

    for (size_t i = 0; i < rez.broj_redova; i++)
        for (size_t j = 0; j < rez.broj_kolona; j++)
            rez.elementi[i][j] = m1.elementi[i][j] + m2.elementi[i][j];

    return rez;

}

// Operator - za oduzimanje matrica
Matrica operator - (const Matrica& m1, const Matrica& m2) {

    // Provjeravamo da li su matrice A i B zapravo skalari u matricnoj formi
    if (m1.skalar && m2.skalar)
        return Matrica(m1[0][0] - m2[0][0]);

    // Provjera podudarnosti matrica
    if (m1.broj_redova != m2.broj_redova || m1.broj_kolona != m2.broj_kolona)
        throw "Nepodudarne dimenzije matrica za oduzimanje!";

    Matrica rez(m1.broj_redova, m1.broj_kolona);

    for (size_t i = 0; i < rez.broj_redova; i++)
        for (size_t j = 0; j < rez.broj_kolona; j++)
            rez.elementi[i][j] = m1.elementi[i][j] - m2.elementi[i][j];

    return rez;

}

// Operator * koji koristi funkciju Strassen za mnozenje dvije matrice
Matrica operator * (const Matrica& A, const Matrica& B) {

    if (A.skalar && B.skalar)
        return Matrica(A[0][0] * B[0][0]);
    else if (A.skalar)
        return A[0][0] * B;
    else if (B.skalar)
        return A * B[0][0];

    // Ako matrice nisu odgovarajuce za mnozenje baca se izuzetak
    if (A.broj_kolona != B.broj_redova)
        throw "Matrice nisu odgovarajuce za mnozenje.";

    // Pozivanje funkcije Strassen
    Matrica C(A.broj_redova, B.broj_kolona);
    C = Strassen(A, B);

    // Uredjivanje rezultata
    C.ispraviNule();

    C.obrisiVisak(A,B);

    return C;
}

Matrica operator * (double skalar, const Matrica& m) {

    // Deklarisemo matricu koju cemo vratit, koja je kopija matrice koja je proslijedjena kao parametar
    Matrica rez(m);

    // Svaki element matrice mnozimo skalarom
    for (size_t i = 0; i < m.broj_redova; i++)
        for (size_t j = 0; j < m.broj_kolona; j++)
            rez[i][j] *= skalar;
    return rez;
}

// Ista funkcija kao i pretthodna napravljena radi komutativnosti
Matrica operator * (const Matrica& m, double skalar) {

    Matrica rez(m);

    for (size_t i = 0; i < m.broj_redova; i++)
        for (size_t j = 0; j < m.broj_kolona; j++)
            rez[i][j] *= skalar;
    return rez;
}

Matrica operator / (const Matrica& A, const Matrica& B) {

    // Provjeravamo da li su matrice A i B zapravo skalari u matricnoj formi
    if (A.skalar && B.skalar)
        return Matrica(A[0][0] / B[0][0]);
    else if (B.skalar)
        return A / B[0][0];
    else
        throw "Matrice se ne mogu dijeliti medjusobno, niti se skalar moze dijeliti matricom.";

}
Matrica operator / (const Matrica& A, double skalar) {

    if (skalar == 0)
        throw "Dijeljenje nulom nije definisano.";

    if (A.skalar)
        return Matrica(A[0][0] / skalar);

    Matrica ret(A);
    for (size_t i = 0; i < ret.broj_redova; i++)
        for (size_t j = 0; j < ret.broj_kolona; j++)
            ret[i][j] /= skalar;
    return ret;


}

//  Operator ^ koji koristi algoritam brzog stepenovanja
Matrica Matrica::operator ^ (int eksponent) const {

    if (eksponent < 0) {
        Matrica temp = Inverzna();
        eksponent *= -1;
        return temp ^ eksponent;
    }
    Matrica rez(broj_redova, broj_kolona);
    Matrica baza = *this;

    while (eksponent > 0) {
        if (eksponent % 2 == 1) {
            rez = rez * baza;
        }
        baza *= baza;
        eksponent /= 2;
    }
    return rez;
}

// Funkcija za brojanje cifri broja
int duzinaBroja(double broj) {


    string s = to_string(broj);
    int k = s.size() - 1;

    while(s[k] == '0') {

        s.erase(s.begin() + k);
        k--;
    }
    return s.size();
}

// Operator  <<  za ispisivanje matrica
ostream& operator  <<  (ostream& ispis, const Matrica& m) {

    // Varijabla koju koristimo da bi nasli broj sa najvise cifara
    int najduzi = 0;
    for (size_t i = 0; i < m.broj_redova; i++)
        for (size_t j = 0; j < m.broj_kolona; j++)
            najduzi = max (najduzi, duzinaBroja(m[i][j]));

    // Ispisujemo matricu postavljajuci sirinu ispisa na duzinu najduzeg broja
    for (size_t i = 0; i < m.broj_redova; i++) {
        for (size_t j = 0; j < m.broj_kolona; j++) {
            ispis  <<  setw (najduzi)  <<  m[i][j] << " ";
        }

        ispis << endl;
    }

    return ispis;
}

// Funkcija koja sluzi za izdvajanje matrice iz ulaznog toka
Matrica uzmiMatricu(istream& unos) {

    // Matrica rezultat
    Matrica m;
    vector<double>red;
    char znak;
    double broj;

    while (1) {

            // Uzimamo znak po znak
            znak = unos.get();
            if (znak == '[')
                znak = unos.get();
            if (isdigit(znak) || znak == '-') {
                unos.putback(znak);
                unos  >>  broj;
                red.push_back(broj);
            }
            else if (znak == ';' || znak == ']') {
                m.elementi.push_back(red);

                // Praznimo red
                red.clear();

                // Kada dodje do ] prestaje sa izdvajanjem matrice
                if (znak == ']')
                    break;
            else if (znak == ' ')
                unos.get();
            }
    }

    // Provjeravamo da li je matrica grbava
    for (size_t i = 0; i < m.elementi.size(); i++)
        if (m.elementi[i].size() != m.elementi[0].size())
            throw "Matrica je grbava";

    // Postavljamo broj redova i kolona matrice
    m.broj_redova = m.elementi.size();
    m.broj_kolona = m.elementi[0].size();

    return m;

}

bool Cifra(char znak) { return isdigit(znak); }

bool Otvorena(char znak) {

    return (znak == '(' || znak == '[' || znak == '{');

}

bool Zatvorena(char znak) {

    return (znak == ')' || znak == ']' || znak == '}');
}

bool Operacija(char znak) {

    return (znak == '+' || znak == '-' || znak == '*' || znak == '/' || znak == '^');
}

Matrica Rezultat(const Matrica& A, const Matrica& B, char operacija) {

    Matrica rez;

    if (operacija == '+')
        rez = A + B;

    else if (operacija == '-')
        rez = A - B;

    else if (operacija == '*')
        rez = A * B;
    else if (operacija == '/')
        rez = A / B;
    else if (operacija == '^') {
        if (!B.getSkalar())
            throw "Nakon operacije '^' mora ici ili 'T' ili skalar.";
        rez = A^B[0][0];
    }

    return rez;

}
char dajOtvorenu(char znak) {

    if (znak == ')')
        return '(';
    else if (znak == ']')
        return '[';
    else if (znak == '}')
        return '{';

    return '0';

}
int Prioritet(char operacija) {

    if (operacija == '+' || operacija == '-')
        return 1;
    else if (operacija == '*' || operacija == '/')
        return 2;
    else if (operacija == '^' || operacija == 'T')
        return 3;
    return 0;

}

//  Funkcija koja prima dva steka te vrsi operaciju
void izvrsiOperaciju(stack<Matrica>&matrice, stack<char>& znakovi) {

    if (znakovi.top() == 'T') {
        Matrica A = matrice.top();
        matrice.pop();
        A = A.Transponovana();
        matrice.push(A);
        znakovi.pop();
        znakovi.pop();
        return;
    }
    if (znakovi.empty() || Prioritet(znakovi.top())==0)
        throw "Izraz nije validan";

    char operacija = znakovi.top();
    znakovi.pop();

    if (matrice.size() < 2)
        throw "Izraz nije validan";

    Matrica B = matrice.top();
    matrice.pop();
    Matrica A = matrice.top();
    matrice.pop();

    Matrica C = Rezultat(A, B, operacija);
    matrice.push(C);

}

istream& operator  >>  (istream& unos, Matrica& m) {

    try{
        stack<Matrica>matrice;
        stack<char>s;

        //  Svi moguci simboli koji se mogu naci u izrazu
        enum Simbol{otvorena, zatvorena, matrica, operacija, T, I};
        Simbol prethodni = otvorena;
        char znak;

        //  Uzimamo prazna mjesta koja se nadju prije izraza
        while (unos.peek() == ' ')
            unos.get();

        //  Sve dok je unos.peek() razlicito od kraja unosa
        while (unos.peek() != '\n') {

            //  Ako je u pitanju orvorena zagrada, dodajemo je na stek s
            if (Otvorena(znak = unos.peek())) {
                if (prethodni == zatvorena || prethodni == matrica)
                    throw "Izmedju zatvorene zagrade i otvorene ili matrice i otvorene mora doci operacija";

                //  Ovdje podrazumijevao da se zagrada '[' koristi samo za unos matrice
                if (znak == '[') {
                    Matrica A = uzmiMatricu(unos);
                    matrice.push(A);
                    prethodni = matrica;
                }

                else {
                    s.push(unos.get());
                    prethodni = otvorena;
                }
            }

            //Ako je u pitanju zatvorena, provjeravamo prethodni te ako je uredu vrsimo operacoje dok ne dodjemo do otvorene
            else if (Zatvorena(unos.peek())) {

                char znak = unos.get();
                if (prethodni == operacija || prethodni == otvorena)
                    throw "Prije zatvorene ne moze biti operacija ili otvorena";

                //  Izvrsavamo operacije
                while ((!s.empty() && Operacija(s.top())) || s.top() == 'T')
                    izvrsiOperaciju(matrice, s);

                if (s.empty() || s.top() != dajOtvorenu(znak))
                    throw "Zagrade nisu dobro balansirane";

                else s.pop();
                prethodni = zatvorena;
            }

            //Ako je u pitanju operacija, provjeravamo prioritet te operacije i one sa vrha steka, ako je prioritet vrha steka veci izvrsavamo operaciju
            else if (Operacija(unos.peek())) {

                if (unos.peek() == '-' && prethodni == otvorena) {
                    double broj;
                    unos  >>  broj;
                    matrice.push(Matrica(broj));
                    prethodni = matrica;
                }

                else if (prethodni == operacija || prethodni == otvorena)
                    throw "Ne moze operacija nakon otvorene zagrade ili druge operacije";

                else{
                    while (!s.empty() && Operacija(s.top()) && Prioritet(unos.peek()) <= Prioritet(s.top()))
                    izvrsiOperaciju(matrice, s);

                    s.push(unos.get());
                    prethodni = operacija;
                }


            }

            //  Ako je u pitanju cifra uzimamo taj broj, te ga pretvarao u skalarnu matricu da bi mogao ici u stek matrica
            else if (Cifra(unos.peek())) {
                if (prethodni == zatvorena)
                    throw "Ne moze broj nakon zatvorene zagrade";

                double broj;
                unos >> broj;
                matrice.push(Matrica(broj));
                prethodni = matrica;
            }

            //  Ako je znak T, dodajemo ga na stek s
            else if (unos.peek() == 'T') {
                if (s.top() != '^')
                    throw "Prije 'T' mora biti '^'.";

                s.push(unos.get());
                prethodni = T;
            }

            // Ako je I, uzimamo znak te broj poslije njega te dodajemo jedinicnu matricu na stek matrica
            else if (unos.peek() == 'I') {
                if (prethodni == zatvorena || prethodni == matrica)
                    throw "Ne mooze prije jedinicne matrica ili zatvorena";

                else {
                    unos.get();
                    int broj;
                    unos  >>  broj;
                    Matrica A(broj, broj);
                    matrice.push(A);
                }

                prethodni = matrica;
            }

            //  Ako je u pitanju razmak samo nastavi dalje
            else if (unos.peek() == ' ') { unos.get(); }

            // Provjerili smo sve moguce znakove, pa ako je neki znak razlicit od tih, on je nepoznat
            else throw "Nepoznat znak";
        }

    //  Izvrasavmo preostale operacije
    while (!s.empty()) {
        izvrsiOperaciju(matrice,s);
    }

    //  Na steku bi trebala ostati samo jedna matrica
    if (matrice.size() != 1)
        throw "Izraz nije dobro napisan";

    // Kopiramo matricu u proslijedjenu
    else m = matrice.top();

    }
    catch(char const poruka[]) {
        cout << poruka;
    }
    catch(...) {}

    return unos;
}
#endif //  MATRICA_CPP

