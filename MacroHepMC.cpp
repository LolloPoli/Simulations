#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include <fstream>
#include <iostream>
#include <cmath>

// MACRO PER CONVERTIRE I DATI HEPMC DI UNA SIMULAZIONE DI SIDIS TRA ELETTRONE E PROTONE 
// CON ANALISI DI UN ADRONE NELLO STATO FINALE ESEGUITA CON PYTHIA6
// _____________________________________________________________________
// per compilare ricordati di inserire la libreria HepMC:  g++ -o Macro6 Macro6.cpp -lHepMC  e poi ./Macro6 output.hepmc
// _____________________________________________________________________
// Nei file HepMC le informazioni sulla simulazione sono divisi in righe, ognuna di esse denominata con diverse lettere
// La riga con E indica l'inizio dell'evento e alcune informazioni di base
// La riga con W indica il peso dell'evento, utile nel caso di correzioni
// La riga con A le informazioni dell'evento (i.e. il numero di particelle)
// La riga con V rappresenta il vertice di interazione e fornisce le sue coordinate
// Le righe con P forniscono: PID, STATUS, M1, M2, Px, Py, Pz, E.
// Lo STATUS sara': 1. particella stabile  2. instabile e puo' decadere  3. generata durante la simulazione
// M1 e M2 indicano gli indici delle particelle madri, indici relativi all'ordine in cui sono elencate le particelle
// _____________________________________________________________________________________________________________________________

// Funzione per calcolare il prodotto scalare tra due quadrivettori
double DotProduct(const HepMC::FourVector& vector1, const HepMC::FourVector& vector2) {
    double dot = vector1.e() * vector2.e() - (vector1.px() * vector2.px() + vector1.py() * vector2.py() + vector1.pz() * vector2.pz());
    return dot;
}

// funzione per Q^2
double QSquared(const HepMC::FourVector& q) {
    double Q2 = -q.m2();
    return Q2;
}
// Calcolo per le variabili di Bjorken
double x_Bjorken(const HepMC::FourVector& q, const HepMC::FourVector& protMom){
    double Q2 = QSquared(q);
    double Pq = DotProduct(protMom, q);
    double x_b = Q2 / (2 * Pq);
    return x_b;
}
double y_Bjorken(const HepMC::FourVector& q, const HepMC::FourVector& eMom_i, const HepMC::FourVector& protMom){
    double Pq = DotProduct(protMom, q);
    double Pk = DotProduct(protMom, eMom_i);
    double y_b = (Pq) / (Pk);
    return y_b;
}
// Funzioni della variabile di SIDIS z_h
double z_hadron(const HepMC::FourVector& q, const HepMC::FourVector& protMom, const HepMC::FourVector& hadMom){
    double PPh = DotProduct(protMom, hadMom);
    double Pq = DotProduct(protMom, q);
    double z_h = (PPh) / (Pq);
    return z_h;
}
// calcolo di W^2
double WSquared(const HepMC::FourVector& q, const HepMC::FourVector& protMom){
    HepMC::FourVector sum = q + protMom;
    double W = sum.m2();
    return W;
}

int main() {

    // RICORDATI DI INSERIRE IL NOME CORRETTO DEL FILE 
    // std::string inputFileName = "output.hepmc";
    std::cout << "Inserisci il nome del file di input (es. output.hepmc): ";
    std::string inputFileName;
    std::cin >> inputFileName;

    // Nome del file in cui vogliamo salvare i nostri dati (.txt o .doc vedi un po' tu)
    std::string outputFileName = "output_analysis.txt";

    // Creazione di un oggetto con IO_AsciiParticles per leggere il file HepMC
    HepMC::IO_AsciiParticles ascii_in(inputFileName, std::ios::in);

    if (ascii_in.failed()) {
        // Gestione dell'errore in caso di fallimento nell'apertura del file
        std::cerr << "Errore: Impossibile aprire il file di input." << std::endl;
        return 1;  
    }

    // Aprire un file di output per la scrittura
    std::ofstream outputFile(outputFileName);

    // Iterare attraverso tutti gli eventi nel file HepMC
    int eventCounter = 0;

    while (!ascii_in.failed()) {
        // Creazione di un oggetto GenEvent per contenere l'evento letto
        //HepMC::GenEvent* evt = ascii_in.read_next_event();
        auto evt = ascii_in.read_next_event();

        if (!evt) {
            // Fine del file
            break;
        }

        // Incrementare il contatore degli eventi
        eventCounter++;

        // Variabile per i momenti, probabilmente non necessari per quello che facciamo
        HepMC::FourVector totalMomentum(0.0, 0.0, 0.0, 0.0);
        // HepMC::FourVector electronMomentum(0.0, 0.0, 0.0, 0.0);
        // HepMC::FourVector protonMomentum(0.0, 0.0, 0.0, 0.0);

        // Considerando che e- sia generato per primo e il protone per secondo
        HepMC::GenParticle* electron = *(evt->particles_begin());
        HepMC::GenParticle* proton = *(evt->particles_begin() + 1);
        // Non sono sicuro serva
        // HepMC::FourVector momentum = particle->momentum();

        // Calcolo dell'energia prima e dopo l'urto dell'elettrone
        HepMC::FourVector eMom_i = electron->momentum();
        double k_i = eMom_i.e();
        HepMC::FourVector eMom_f= (electron->end_vertex()->particles_out_const_begin())->momentum();
        double k_f = eMom_f.e();
        // HepMC::FourVector q = (k_i - k_f, eMom_i.px() - eMom_f.px(), eMom_i.py() - eMom_f.py(), eMom_i.pz() - eMom_f.pz());
        // piu' compatto
        HepMC::FourVector q = eMom_i - eMom_f;
        // NEL CASO USASSI QUESTA DEVO MODIFICARE LE FUNZIONI CHE SONO TUTTE SU k_i e k_f, E METTERLE SU eMom_i e eMom_f

        // prendo Q2 dalla funzione
        double Q2 = QSquared(q);

        // Calcolo del quadrimomento del protone, passggio non necessario ma mettiamolo per comprensione
        HepMC::FourVector protMom = proton->momentum();
        // uso .e() per farmi restituire l'energia dal quadrivettore 
        double E_p = protMom.e();
        double Px = protMom.px();
        double Py = protMom.py();
        double Pz = protMom.pz();
        double P_tot = std::sqrt(Px*Px + Py*Py + Pz*Pz);

        // Calcolo delle variabili di Bjorken tramite le funzioni
        double x_B = x_Bjorken(q, protMom);
        double y_B = y_Bjorken(q, protMom, eMom_i);

        // Prendiamo anche la variabile s e W2
        double s = std::pow((k_i + E_p), 2);
        double W2 = WSquared(q, protMom);

        // definisco la beta per il boost in Breit frame
        double beta = - (q + protMom) / (k_i + E_p);

        outputFile<< "__________________________________________________________________________________________" << std::endl;
        outputFile << " " << std::endl;
        outputFile << "Event #" << eventCounter << std::endl;
        outputFile<< "Q^2: " << Q2 << std::endl;
        outputFile << "Bjorken variables: x_B = " << x_B << "  " << "y_B = " << y_B << std::endl;
        outputFile << "The invariant mass W of the final state is: " << W2 << std::endl;


        // Ciclo su tutte le particelle generate
        for (HepMC::GenEvent::particle_iterator it = evt->particles_begin(); it != evt->particles_end(); ++it) {
            HepMC::GenParticle* particle = (*it);

            /* Controlla se la particella è stata generata durante l'interazione
            if (particle->production_vertex() && particle->production_vertex()->particles_in_const_begin() != particle->production_vertex()->particles_in_const_end()) {
                // Ottenere informazioni sulla particella generata
                int pdgId = particle->pdg_id();
                HepMC::FourVector momentumGenerated = particle->momentum();   
                outputFile << "Particella generata identificabile con PDG ID: " << pdgId << std::endl;
            } */

            // Ottenere informazioni sulla particella generata
            int pdgId = particle->pdg_id();
            HepMC::FourVector Momentum = particle->momentum(); 
            // HepMC::FourVector MomentumBoosted = Momentum;  // per il boost 

            // VOLESSI EFFETTUARE UN BOOST NEL SISTEMA DI BREIT (cosi metto l'energia del fotone a zero eheh)
            // MomentumBoosted.boost(0.0, 0.0, q); // boost lungo z
            // MomentumBoosted.boost(q*beta, 0.0, 0.0); // boost lungo E


            // Aggiungere il momento della particella al momento totale
            // totalMomentum += momentum;  
            // pero' potrebbe non essere sensato siccome siamo in SIDIS quindi non consideriamo effettivamente il sistema completo...

            // Qui dipende che particella vogliamo analizzare, proviamo con un pione carico 
            if (pdgId == 211 || pdgId == -211) {
                outputFile << "E' stato generato un pione (" << pdgId << ")"<< std::endl;
                // estapoliamo il momento del nostro adrone
                HepMC::FourVector hadMom = particle->momentum();
                double E_h = hadMom.e();

                // ora possiamo calcolare la nostra z_h
                double z_h = z_hadron(q, protMom, hadMom); 
                // e' importante siccome per le FF e' analoga a x_B

                outputFile << "P_h = " << protMom << std::endl;
                outputFile << "z_h = " << z_h << std::endl;
            }

            // potrebbe servirci saperlo nel caso nessun pione venisse generato? lo scopriremo
            else {
                outputFile << "Nessun pione e' stato generato in questa collisione" << std::endl;
            } 
            // per avere una maggior chiarezza nella stampa
            outputFile << " " << std::endl;
            outputFile << "____________________________________________________________________________________" << std::endl;

        }
        // Elimina l'evento a fine ciclo
        delete evt;
    }
    // Chiudo il file di output
    outputFile.close();

    return 0;
}