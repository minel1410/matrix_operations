#ifndef MATRICA_H
#define MATRICA_H
#include<iostream>
#include<vector>
#include<cmath>
#include<stack>
#include<iomanip>
using namespace std;

class Matrica{

    size_t broj_redova, broj_kolona;                            // S obzirom da ove dvije varijable predstavljaju velicinu matrice, korisitimo size_t
    vector<vector<double>>elementi;                             // Koristimo vektor vektora da cuvamo elemente matrice
    bool skalar = false;                                        // Svaka matrica sadrzi bool, koji nam govori je li to prava matrica ili je samo skalar u matricnoj formi

    bool jeLiKvadratna() const ;                                // Provjera da li je matrica kvadratna
    Matrica napraviKvadratnu () const;                          // Funkcija koja uzima matricu, te vraca kvadratnu matricu tako sto kopira originalnu te dodaje nula redove i kolone
    Matrica napraviStepenOdDva () const;                        // Funkcija koja pravi kvadratnu, ali takvu da je broj redova stepen od 2

    void obrisiVisak (const Matrica&, const Matrica&);          // Funkcija koja mice nula redove i kolone koji su nastali u Strassenovom algoritmu
    friend Matrica Strassen (const Matrica&, const Matrica&);   // Funkcija koja mnozi dvije matrice u O(n^2.81)
    void ispraviNule ();                                        // Funkcija koja zaokruzuje jako male brojeve na 0, koji nastaju zbog floating point aritmetike

public:

    Matrica ();                                                 // Konstruktor bez parametara
    Matrica (int, int);                                         // Jedinicna matrica dimenzija m x n
    Matrica (const double&);                                    // Konstruktor koji predstavlja skalar u matricnoj formi, postavlja skalar = true
    Matrica& operator = (const double&);                        // Operator dodjele za koji prima skalar
    Matrica (const Matrica&, int, int);                         // Konstruktor koji kopira matricu ali bez reda i kolone koji se prosljedjuju
    Matrica (const Matrica&);                                   // Konstruktor kopije
    Matrica (Matrica&&);                                        // Pomjerajuci konstruktor kopije
    Matrica& operator = (const Matrica&);                       // Operator dodjele
    Matrica& operator = (Matrica&&);                            // Pomjerajuci operator dodjele
    ~Matrica ();                                                // Destruktor

    Matrica Transponovana () const;                             // Konstantna fja koja vraca transponovanu matricu u odnosu na prosljedjenju
    Matrica Inverzna () const;                                  // Konstantna fja koja vracu inverznu matricu

    bool getSkalar () const;                                    // Geter za skalar

    void dodajRedNaKraj (vector<double>);                       // Funkcija koja dodaje prosljedjeni vektor kao red na kraj matrice
    void dodajRedNaPocetak (vector<double>);                    // Funkcija koja dodaje prosljedjeni vektor kao red na pocetak matrice
    void dodajKolonuNaKraj (vector<double>);                    // Funkcija koja prosljedjeni vektor dodaje kao kolonu na kraj matrice
    void dodajKolonuNaPocetak (vector<double>);                 // Funkcija koja prosljedjeni vektor dodaje kao kolonu na pocetak matrice

    double Determinanta () const ;                              // Racunanje determinante

    vector<double>& operator [] (int);                          // Operator koji vraca vektor iz vektora elementi po referenci, gdje je matrica nekonstantna
    const vector<double>& operator [] (int) const;              // Isti operator koji vraca konstantnu vrijednost iz konstantne matrice

    bool operator == (const Matrica&);                          // Poredbeni operatrii
    bool operator != (const Matrica&);

    void operator += (const Matrica&);
    void operator -= (const Matrica&);
    void operator *= (const Matrica&);
    void operator *= (double);
    void operator /= (double);
    Matrica operator ^ (int) const;                             // Brzo stepenovanje

    friend Matrica operator + (const Matrica&, const Matrica&); // Operator +
    friend Matrica operator - (const Matrica&, const Matrica&); // Operator -
    friend Matrica operator * (const Matrica&, const Matrica&); // Operator * koji koristi Strassenov algoritam
    friend Matrica operator * (double, const Matrica&);         // Mnozenje matrice skalarom
    friend Matrica operator * (const Matrica&, double);
    friend Matrica operator / (const Matrica&, const Matrica&); // Dijeljenje matrica, ali samo kada je druga ili obje matrice skalar
    friend Matrica operator / (const Matrica&, double);         // Dijeljenje matrica, ali samo kada je druga ili obje matrice skalar

    friend ostream& operator << (ostream&, const Matrica&);     // Operator za ispis matrice
    friend istream& operator >> (istream&, Matrica&);           // Upis matrice iz ulaznog toka
    friend Matrica uzmiMatricu (istream&);                      // Funkcija koja vraca matricu izdvajanjem iz ulaznog toka
};

#endif //  MATRICA_H
