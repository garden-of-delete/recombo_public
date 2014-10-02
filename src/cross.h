/* fichier croisements_l.h
 * date creation	:	25/06/97
 *	
 * description		:
 *	Fichier entete a utiliser avec le module croisements_l.c
 *
 *	Auteur		MEUNIER-GUTTIN-CLUZEL Siegfried
 *				LESP-UMR 6614 Coria- Operation 17
 *
 *	historique		:
 *	- issu de gaussred1.c
 *	- issu de croisements1.h
 *	- (17/06/97)	ajout des codes de gauss signes
 *	- (25/06/97)	version pour les liens
 *	- (9/10/97)		mise en forme
 *	- (10/06/98)	nouvelle gestion des erreurs
 *	- (6/07/98)		ajout fonctions
 *							suivant_transverse_de(),
 *							lit_croisement(),lit_croisement2()
 *							saute_triviaux()
 *
 *	version du 6/07/98
 **************************************************************************/
#ifndef		CROISEMENTS


#define		CROISEMENTS
#define		DIRECT			0
#define		INVERSE			1
#define		DEBUT			0
#define		MOINS			-1
#define		PLUS			1
#define		ABOVE			0
#define		BELOW			1
#define		NON				0
#define		OUI				1
#define		MAX_NUM			50
#define		MAX_SEQUENCE	MAX_NUM*2
#define		MAX_MEMBRES		10
#define		ECRAN			0
#define		LATEX			1
#define		MAPPLE			2
#define		BINAIRE			3

#define		ERREUR_SYNTAXE	-1
#define		ERREUR_MEMOIRE	-2
#define		ERREUR_FIN		-3

#define ERROR_RGS_1   -101
#define ERROR_RGS_2   -102
#define ERROR_RGS_3   -103
#define ERROR_RGS_4   -104
#define ERROR_RGS_5   -105
#define ERROR_RGS_6   -106
#define ERROR_RGS_7   -107
#define ERROR_RGS_8   -108
#define ERROR_RGS_9   -109
#define ERROR_RGS_10  -110
#define ERROR_RGS_11  -111
#define ERROR_RGS_12  -112
#define ERROR_RGS_13  -113
#define ERROR_RGS_14  -114


#include	<stdlib.h>
#include	<stdio.h>

typedef struct CR {
	int			numero;
	char		sens;
	struct CR	*precedent,
				*suivant;
#ifdef	CODES_SIGNES
	char		signe;
#endif
}CROISEMENT;
extern	char	option_affichage;
			
/*****************************************************************************
 * FONCTIONS DIVERSES DE GESTION DES CHAINES DE CROISEMENTS
 *****************************************************************************/
/*** recherche le croisement suivant en sautant eventuellement le debut *****/
CROISEMENT
*suivant_de();
/*** recherche le croisement precedent en sautant eventuellement le debut *****/
CROISEMENT
*predecesseur_de();
/*** saute si c'est le croisement factice de debut ***************************/
CROISEMENT
*saute_debut();
/*** renvoie OUI si le croisement p est le croisement factice de debut *******/
char
est_debut();
/*** retourne le numero ( entre 0 et maxNum ) d'un croisement ***/
int
numero();
#define		NUMERO( pc )	((pc)->numero)
/*** retourne le prefixe ( ABOVE ou BELOW ) d'un croisement *******/	
char
prefixe();
#define		PREFIXE( pc )	((pc)->sens)
/*** retourne le signe ( PLUS ou MOINS ) d'un croisement *******/
#ifdef	CODES_SIGNES	
char
signeCroisement();
#define		SIGNE_CROISEMENT( pc )	((pc)->signe)
#endif
/***************** croisement suivant en prenant l'arc transverse *****/
CROISEMENT
*suivant_transverse_de();
/***************** lit les parametres d'un croisement ********************/
int
lit_croisement();
int
lit_croisement2();
/***************** saute les noeuds triviaux *****************************/
CROISEMENT
*saute_triviaux();
/************ trouve l'autre croisement de meme numero *******************/
CROISEMENT
*trouve_l_autre();
/*** creation dynamique d'un nouveau croisement *******************/
CROISEMENT
*creerCroisement();
/** etablit les connexions entre 2 croisement consecutifs ******************/
void
connecteCroisements();
/*********** ajoute un croisement a une liste dont on indique le dernier ***/
CROISEMENT
*ajouteCroisement();
/*********** insere un croisement dans une liste apres l'element pointe ***/
CROISEMENT
*insereCroisement();
/************ determine le nombre de croisements (meme fictifs) d'un code*/
int
tailleDuCode();
/************ supprime un croisement d'une sequence *********************/
void
libereCroisement();
/************ supprime completement une sequence de croisements *********/
void
libereCode();
/*********** affiche un croisement sous la forme sn ( s=a ou s=b ) *********/
void
afficheCroisement();

/*** signalement d'une erreur de syntaxe sur la chaine de depart ***********/
void
erreurSyntaxe();
/*** affichage simple du code **( exemple : Ea1b2a3b1a2b3 ) ****************/
void
afficheCode();
/*** entree manuelle du code sous forme d'une chaine de caracteres ****/	
/* attention changement parametres (11/06/98)                         */	
int
entrerCode();
int
entrerCodeDynamique();


#endif
/** fin du ifndef CROISEMENTS **/

FILE *p;
