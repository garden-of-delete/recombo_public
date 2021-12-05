/* fichier polynomeAZ.h
 * description		:
 *	Fichier entete a utiliser avec le module polynomeAZ.c
 *
 *	Auteur		MEUNIER-GUTTIN-CLUZEL Siegfried
 *				LESP-UMR 6614 Coria- Operation 17
 *
 *********************************************************************/

#ifndef		POLYNOMES_AZ

#define		POLYNOMES_AZ

#define			MAX_POLY	50

typedef struct{
	float	coeff;
	int		degA, degZ;
}MONOME;

int polyCopie(MONOME p1[],MONOME p2[], int n2);
int polyProduit(MONOME p1[], int n1,MONOME p2[], int n2,MONOME p3[], int n3);
int polySomme(MONOME p1[], int n1,MONOME p2[], int n2,MONOME p3[], int n3);
int polySimplifie(MONOME p1[], int n1);
void polyAffiche(MONOME p1[], int n1, FILE *pout);
void polyLatex(MONOME p1[], int n1);
void polyMapple(MONOME p1[], int n1);

#endif
