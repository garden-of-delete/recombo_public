/* fichiers croisements_l.c
* date creation	:	25/06/97
*	
* description		:
*	Module contenant les fonctions de manipulation de liens pour les
*	programmes tournant autour de la theorie des noeuds. Utilisation
*	des codes de Gauss ( signes ou non suivant que la constante CODES_SIGNES
*	est definie ou non.
*

Module containing functions for manipulating links
programs revolving around the theory of knots. use
codes of Gauss (or not following signs that the constant CODES_SIGNES
is defined or not.

*	Auteur		MEUNIER-GUTTIN-CLUZEL Siegfried
*				LESP-UMR 6614 Coria- Operation 17
*
*	historique		:
*	- issu de gaussred1.c
*	- issu de croisements.c
*		differences au niveau des fonctions suivantes :
*			trouve_l_autre(p1)	->	trouve_l_autre(p1,lien)
*			tailleDuCode(code)	->	tailleDuCode(lien)
*			libereCode(code)	->	libereCode(lien,nm)
*			afficheCode(code)	->	afficheCode(lien)
*			entrerCode(code,bavard)->	entrerCode(lien,bavard)
*			entrerCodeDynamique(code,chaine,maxChaine,maxTaille,bavard) ->
*			entrerCodeDynamique(lien,chaine,maxChaine,maxTaille,maxMembres,bavard)
*		(6/07/98)
*			suivant_transverse_de(p) -> suivant_transverse_de(p,lien)
*			saute_triviaux(lien,pNumMembre)
*		
*	- (17/06/97)	ajout des codes de gauss signes
*	- (25/06/97)	gestion des liens
*	- (30/06/97)	verification des codes
*	- (03/06/98)	ajout entree dynamique ( pas de limite de taille )
*	- (16/06/98)	correction entrerCode pour gerer le noeuds triviaux
*	- (6/07/98)		ajout fonctions
*						suivant_transverse_de(),saute_triviaux()
*							lit_croisement(),lit_croisement2()
*
*	Compilation	cc -c croisements_l.c
*	exemples 	b1-a2-|a3+b4+b2-a1-|a4+b3+				(avec signes)
*				b1a2|a3b4b2a1|a4b3						(sans signes)
*
*	version du (6/07/98)
* History:
* - Derived from gaussred1.c
* - Derived from croisements.c
* Differences in the following functions:
Trouve_l_autre * (p1) -> trouve_l_autre (p1, link)
* TailleDuCode (code) -> tailleDuCode (link)
* LibereCode (code) -> libereCode (link, nm)
* AfficheCode (code) -> afficheCode (link)
* EntrerCode (code, verbose) -> entrerCode (link, talkative)
* EntrerCodeDynamique (code, string, maxChaine, maxTaille, talkative) ->
* EntrerCodeDynamique (link, chain, maxChaine, maxTaille, maxMembres, talkative)
* (7.6.98)
Suivant_transverse_de * (p) -> suivant_transverse_de (p, link)
* Saute_triviaux (link pNumMembre)
*
* - (06/17/97) adding signs codes gauss
* - (25/06/97) link management
* - (06/30/97) verification codes
* - (03/06/98) Adding dynamic entry (no size limit)
* - (06/16/98) Correction entrerCode to manage the trivial knots
* - (07/06/98) Added functions
Suivant_transverse_de * (), saute_triviaux ()
Lit_croisement * (), lit_croisement2 ()
*
* Compilation cc-c croisements_l.c
* Examples b1-a2-| a3 + b4 + b2-a1-| a4 + b3 + (with signs)
* B1a2 | a3b4b2a1 | a4b3 (without signs)
*
* Version (6.7.98)
**************************************************************************/
#define			CODES_SIGNES
#include		"cross.h"
#include    <string.h>
char			option_affichage;
/*****************************************************************************
* FONCTIONS DIVERSES DE GESTION DES CHAINES DE CROISEMENTS POUR LIENS
*****************************************************************************/
/************************************************* ****************************
* OTHER FUNCTIONS OF MANAGEMENT FOR CHAIN ​​LINKS CROSSINGS
************************************************** ***************************/
extern void report_error (char *, int, char *);
extern void report_error2 (char *, int, char *, int);

/*** recherche le croisement suivant en sautant eventuellement le debut *****/
/*** Search the next crossing possibly jumping the start *****/
CROISEMENT
*suivant_de(p)
CROISEMENT		*p;
{
  p = p->suivant;
  if(p->numero == DEBUT)
    p = p->suivant;
  return p;
}

/*** recherche le croisement precedent en sautant eventuellement le debut *****/
/*** Previous research crossing possibly jumping the start *****/
CROISEMENT
*predecesseur_de(p)
CROISEMENT		*p;
{
  p = p->precedent;
  if(p->numero == DEBUT)
    p = p->precedent;
  return p;
}

/*** saute si c'est le croisement factice de debut ***************************/
/*** Skip if this is the beginning of artificial crossing ***************************/
CROISEMENT
*saute_debut(p)
CROISEMENT		*p;
{
  if(p->numero == DEBUT)
    p = p->suivant;
  return p;
}

/*** renvoie OUI si le croisement p est le croisement factice de debut *******/
/*** Returns YES if the crossing is the crossing p dummy start *******/
char
est_debut(p)
CROISEMENT		*p;
{
  return (p->numero == DEBUT ? OUI : NON );

}

/*** retourne le numero ( entre 0 et maxNum ) d'un croisement ***/
int
numero(pc)
CROISEMENT		*pc;
{
  return pc->numero -1;
}
/*** retourne le prefixe ( ABOVE ou BELOW ) d'un croisement *******/	
char
prefixe(pc)
CROISEMENT		*pc;
{
  return pc->sens;
}
/*** retourne le signe ( PLUS ou MOINS ) d'un croisement *******/
#ifdef	CODES_SIGNES	
char
signeCroisement(pc)
CROISEMENT		*pc;
{
  return pc->signe;
}
#define		SIGNE_CROISEMENT( pc )	((pc)->signe)
#endif

/***************** croisement suivant en prenant l'arc transverse *****/
/* version liens (26/06/97)											  */
CROISEMENT
*suivant_transverse_de(p,lien)
CROISEMENT		*p, **lien;
{
  CROISEMENT		*p1;

  p1= trouve_l_autre(p,lien);
  return suivant_de(p1);
}
/***************** lit les parametres d'un croisement ********************/
int
lit_croisement(p,pPrefixe,pNumero,pSigne)
CROISEMENT		*p;
char			*pPrefixe,*pSigne;
int				*pNumero;
{
  if(p==NULL)return ERROR_RGS_8;
  *pPrefixe= prefixe(p);
  *pNumero= numero(p);
  *pSigne= signeCroisement(p);
  return *pNumero*2 + *pPrefixe;
}
/***************** lit les parametres d'un croisement ********************
* retourne un numero reduit en utilisant un tableau 'dictionnaire'
*	numTab, tel que numeroReduit = numTab[numero(p)]
*************************************************************************/
int
lit_croisement2(p,pPrefixe,pNumeroReduit,pSigne,numTab)
CROISEMENT		*p;
char			*pPrefixe,*pSigne;
int				*pNumeroReduit;
short 			*numTab;
{
  if(p==NULL || numTab==NULL)return ERROR_RGS_9;
  *pPrefixe= prefixe(p);
  *pNumeroReduit= numTab[numero(p)];
  *pSigne= signeCroisement(p);
  return *pNumeroReduit*2 + *pPrefixe;
}

/***************** saute les noeuds triviaux ****************************
*	mod. 1/07/98
*	retourne NULL et numMembre = ERREUR_MEMOIRE si il y a un probleme de
*	pointeur, si numMembre >= 0 alors on n'a pas trouve de noeuds non trivial
************************************************************************/
CROISEMENT
*saute_triviaux(lien,pNumMembre)
CROISEMENT		**lien;
int				*pNumMembre;
{
  int				nm;
  CROISEMENT		*pDeb;

  if( lien==NULL || lien[0]==NULL ){
    *pNumMembre = ERROR_RGS_10;
    return NULL;
  }
  nm = *pNumMembre;
  if( lien[nm]==NULL ) return NULL;

  pDeb=saute_debut(lien[nm]);
  while( est_debut(pDeb) && lien[nm+1]!=NULL ){
    nm++;
    pDeb=saute_debut(lien[nm]);
  }

  *pNumMembre = nm;
  if(est_debut(pDeb))return NULL;		/* pas trouve de noeud non trivial */
  return pDeb;
}

/*** creation dynamique d'un nouveau croisement *******************/
/*** version sans signes ******************************************/
#ifndef	CODES_SIGNES
CROISEMENT
*creerCroisement(numero,sens)
int			numero;
char		sens;
{
  CROISEMENT	*ptC;

  ptC = (CROISEMENT *)malloc(sizeof(CROISEMENT));
  if(ptC != NULL ){
    ptC->numero = numero;
    ptC->sens = sens;
    ptC->precedent = NULL;
    ptC->suivant = NULL;
  }
  return (ptC);
}
#endif
/*** version avec signes ******************************************/
#ifdef	CODES_SIGNES
CROISEMENT
*creerCroisement(numero,sens,sig)
int			numero;
char		sens,sig;
{
  CROISEMENT	*ptC;

  ptC = (CROISEMENT *)malloc(sizeof(CROISEMENT));
  if(ptC != NULL ){
    ptC->numero = numero;
    ptC->sens = sens;
    ptC->signe = sig;
    ptC->precedent = NULL;
    ptC->suivant = NULL;
  }
  return (ptC);
}
#endif
/** etablit les connexions entre 2 croisement consecutifs ******************/
void
connecteCroisements(p1,p2)
CROISEMENT	*p1, *p2;
{
  if(p1==NULL || p2==NULL ) return;
  p1->suivant = p2;
  p2->precedent = p1;
}

/*********** ajoute un croisement a une liste dont on indique le dernier ***/
/*** version sans signes ******************************************/
#ifndef	CODES_SIGNES
CROISEMENT
*ajouteCroisement(dernier,numero,sens)
CROISEMENT	*dernier;
{
  CROISEMENT	*nouveau;		

  if((nouveau = creerCroisement(numero,sens))==NULL){
    report_error (__FILE__, __LINE__,"ERREUR : allocation memoire\n");
    exit(-1);
  }
  connecteCroisements(dernier,nouveau);
  return nouveau;
}
#endif
/*** version avec signes ******************************************/
#ifdef	CODES_SIGNES
CROISEMENT
*ajouteCroisement(dernier,numero,sens,sig)
CROISEMENT	*dernier;
{
  CROISEMENT	*nouveau;		

  if((nouveau = creerCroisement(numero,sens,sig))==NULL){
    report_error (__FILE__, __LINE__,"ERREUR : allocation memoire\n");
    exit(-1);
  }
  connecteCroisements(dernier,nouveau);
  return nouveau;
}
#endif

/*********** insere un croisement dans une liste apres l'element pointe ***/
/*** version sans signes ******************************************/
#ifndef	CODES_SIGNES
CROISEMENT
*insereCroisement(dernier,numero,sens)
CROISEMENT	*dernier;
{
  CROISEMENT	*nouveau;		

  if((nouveau = creerCroisement(numero,sens))==NULL){
    report_error (__FILE__, __LINE__,"ERREUR : allocation memoire\n");
    exit(-1);
  }
  connecteCroisements(nouveau,dernier->suivant);
  connecteCroisements(dernier,nouveau);
  return nouveau;
}
#endif
/*** version avec signes ******************************************/
#ifdef	CODES_SIGNES
CROISEMENT
*insereCroisement(dernier,numero,sens,sig)
CROISEMENT	*dernier;
{
  CROISEMENT	*nouveau;		

  if((nouveau = creerCroisement(numero,sens,sig))==NULL){
    report_error (__FILE__, __LINE__,"ERREUR : allocation memoire\n");
    exit(-1);
  }
  connecteCroisements(nouveau,dernier->suivant);
  connecteCroisements(dernier,nouveau);
  return nouveau;
}
#endif
/************ trouve l'autre croisement de meme numero *******************
*				version pour liens 		(mod.25/06/97)					 */
CROISEMENT
*trouve_l_autre(p1,lien)
CROISEMENT		*p1,**lien;
{
  CROISEMENT		*p2,*pDeb;
  int				n1,nm;
  char			trouve=NON;

  if(p1==NULL || lien==NULL) return NULL;

  n1 = numero(p1);
  p2= suivant_de(p1);
  nm=0;
  while(lien[nm]!=NULL && !trouve){
    pDeb=saute_debut(lien[nm]);
    p2=pDeb;
    do{
      if(numero(p2)==n1 && p2!= p1) trouve=OUI;
      else p2 = suivant_de(p2);
    }while(!trouve && p2!=pDeb);
    if(!trouve){
      nm++;
    }
  }	
  if( trouve ) return p2;
  return NULL;
}


/************ determine le nombre de croisements (non fictifs) d'un code	*
*				version pour liens 		(mod.25/06/97)						*
*			- on ne tient pas compte des croisement fictifs					*/
int
tailleDuCode(lien)
CROISEMENT		**lien;
{
  CROISEMENT		*p;
  int				i=0,nm=0;

  if(lien==NULL)return 0;
  while(lien[nm]!=NULL){
    p=saute_debut(lien[nm]);
    do{
      if(!est_debut(p)) i++;
      p=p->suivant;
    }while(!est_debut(p));
    nm++;
  }

  return (i);
}
/************ supprime un croisement d'une sequence *********************/
void
libereCroisement(p)
CROISEMENT		*p;
{
  if( p==NULL || (p->numero)==DEBUT ) return;
  connecteCroisements(p->precedent,p->suivant);
  free(p);

}
/************ supprime completement une sequence de croisements **********
*				version pour liens 		(mod.25/06/97)					 *
*		- supprime le membre nm du lien									 *
*		- il y a reorganisation du tableau lien							 */
void
libereCode(lien,nm)
CROISEMENT		**lien;
int				nm;
{
  CROISEMENT		*p;

  if(lien==NULL || lien[nm]==NULL) return;
  p=saute_debut(lien[nm]);
  while(!est_debut(p)){
    libereCroisement(p);
    p=p->suivant;
  }
  if(est_debut(p))free(p);
  while(lien[nm+1]!=NULL){
    lien[nm]=lien[nm+1];
    nm++;
  }
  lien[nm]=NULL;
}


/*********** affiche un croisement sous la forme sn ( s=a ou s=b ) *********/
void
afficheCroisement(ptf,p)
FILE		*ptf;
CROISEMENT	*p;
{
  if(est_debut(p))fputc('E',ptf);
  else{
    fputc((p->sens==ABOVE)?'a':'b',ptf);
    fprintf(ptf,"%d",p->numero);
#ifdef	CODES_SIGNES
    fputc((p->signe==PLUS)?'+':'-',ptf);
#endif
  }
}
/*** signalement d'une erreur de syntaxe sur la chaine de depart ***********/
void
erreurSyntaxe(chaine,i)
char		*chaine;
int			i;
{
  int			j;

  for(j=0;j<=i;j++)fputc(chaine[j],stderr);
  fputc('\n',stderr);
  if(i<80){
    for(j=0;j<i;j++)fputc(' ',stderr);
    fputc('^',stderr); fputc('\n',stderr);
  }

}
/*** affichage simple du code ********************************************
*				version pour liens 		(mod.25/06/97)					 *
*	exemples:	b1-a2-|a3+b4+b2-a1-|a4+b3+			(avec signes)		 *
*				b1a2|a3b4b2a1|a4b3					(sans signes)		 */
void
afficheCode(lien)
CROISEMENT	**lien;
{
  CROISEMENT	*p;
  int			nm=0;

  if(lien==NULL) return;
  if(option_affichage == LATEX) putchar('$');

  while((p=lien[nm])!=NULL){
    if(est_debut(p)){
      if(nm>0)putchar('|');
      p = p->suivant;
    }
    while(!est_debut(p)){
      afficheCroisement(stdout,p);
      p = p->suivant;
    }
    nm++;
  }
  if(option_affichage == LATEX){
    putchar('$'); putchar('\n');
  }
  putchar('\n');

}
/*** entree manuelle du code sous forme d'une chaine de caracteres *******
*				version pour liens 		(mod.25/06/97)					 
*	(12/06/98)	ajout parametre bavard + modif gestion erreurs
*	(16/06/98)	correction pour gerer le noeuds triviaux
***************************************************************************/
int
entrerCode(lien,bavard)
CROISEMENT		**lien;
char			bavard;
{
  char		chaine[256],buff[8],nb[MAX_NUM],si[MAX_NUM],pr[MAX_NUM],
    sens,sig,debutMembre=OUI,ttest=NON;
  int			l,i,j,k,nbC,nm,numero;
  CROISEMENT	*debut, *dernier;


  if(lien==NULL)return ERROR_RGS_11;
  for(nm=0;nm<MAX_MEMBRES;nm++)lien[nm]=NULL;
  for(j=0;j<MAX_NUM;j++) nb[j]=0;	


  if(scanf("%s",chaine)==EOF)return  ERREUR_FIN;
  l = strlen(chaine);

  /************** traitement *********************************************/
  nm = 0; 											/* numero de membre */
  i = 0;											 /* numero de caractere */
  nbC = 0;										/* numero de croisement */

  /* saute les caracteres de d\E9but *******************/
  while((chaine[i]=='\t' || chaine[i]==' ' 
    || chaine[i]=='e' || chaine[i]=='E')&&(i<l) )i++;

  /**** boucle principale de lecture et de decodage ************************/
  do{
    if(debutMembre){							/* creation debut membre */
      if(ttest)printf("Allocation membre %d\n",nm);
      if((debut = creerCroisement(DEBUT,-1))==NULL){
        if(bavard)report_error (__FILE__, __LINE__,"ERREUR : allocation memoire\n");
        return ERROR_RGS_12;
      }
      lien[nm]=debut;
      dernier = debut;	
      debutMembre=NON;
    }

    if(chaine[i]=='|' ){ 						 /** debut d'un membre **/
      debutMembre=OUI;
      connecteCroisements(dernier,debut);				   /* on boucle */

      nm++;
      if(nm>=MAX_MEMBRES){
        if(bavard)report_error (__FILE__, __LINE__,"ERREUR : trop de membres\n");
        return ERROR_RGS_13;
      }
      i++;
    }

    /*** decode un croisement **********************************************/
    else {
      if(chaine[i]=='a') sens = ABOVE;						/* prefixe */
      else  if(chaine[i]=='b') sens = BELOW;
      else {
        if(bavard){
          report_error (__FILE__, __LINE__,"ERREUR : caractere illegal :\n");
          erreurSyntaxe(chaine,i);
        }
        return ERREUR_SYNTAXE;
      }

      if(chaine[i+1]<'0' || chaine[i+1]>'9'){					 /* numero */
        if(bavard){
          report_error (__FILE__, __LINE__,"ERREUR : caractere illegal :\n");
          erreurSyntaxe(chaine,i+1);
        }
        return ERREUR_SYNTAXE;
      }
      j=i+1; k=0;
      while( chaine[j]>='0' && chaine[j]<='9' ) buff[k++]=chaine[j++];
      buff[k]='\0';
      numero=atoi(buff);
      i=j;

#ifdef	CODES_SIGNES
      if(chaine[i]=='+') sig = PLUS;							  /* signe */
      else  if(chaine[i]=='-') sig = MOINS;
      else {
        if(bavard){
          report_error (__FILE__, __LINE__,"ERREUR : caractere illegal :\n");
          erreurSyntaxe(chaine,i);
        }
        return ERREUR_SYNTAXE;
      }
      i++;
#endif

      /**** verifications ***************************************/
      nb[numero]++;
      if(nb[numero]>2){
        if(bavard){
          report_error2 (__FILE__, __LINE__,"ERREUR : Trop de passage sur %d\n",numero);
          erreurSyntaxe(chaine,i-2);
        }
        return ERREUR_SYNTAXE;
      }
      if(nb[numero]==1){						/* 1ere visite   */
        si[numero]=sig;
        pr[numero]=sens;
      }
      if(nb[numero]==2){						/* 2eme visite  */
        if( si[numero]!=sig ){
          if(bavard){				
            report_error (__FILE__, __LINE__,"ERREUR : Signe incoherent \n");
            erreurSyntaxe(chaine,i-1);
          }
          return ERREUR_SYNTAXE;
        }
        if( pr[numero]== sens){
          if(bavard){				
            report_error2 (__FILE__, __LINE__,
              "ERREUR : Deux fois meme sens sur %d\n",numero);
            erreurSyntaxe(chaine,i-2);
          }
          return ERREUR_SYNTAXE;
        }
      }
      /**** ajout du croisement dans la list chainee **********************/
      nbC++;
#ifndef	CODES_SIGNES
      if(ttest)printf("Croisement %d sens %c numero %d\n",nbC+1,
        (sens==ABOVE)?'a':'b', numero); 
      dernier = ajouteCroisement(dernier,numero,sens);
#else
      if(ttest)printf("Croisement %d sens %c numero %d signe %c\n",nbC+1,
        (sens==ABOVE)?'a':'b', numero,(sig==PLUS?'+':'-')); 
      dernier = ajouteCroisement(dernier,numero,sens,sig);
#endif

    }									/** fin boucle decodage croisement **/
  }while(i<l);									 /** fin boucle de lecture **/


  connecteCroisements(dernier,debut);						/* on boucle */
  return nbC;
}

/*** entree manuelle du code sous forme d'une chaine de caracteres ****/
/* version dynamique ( pas de limite de taille, version liens )
*	chaine \E0 fournir de longueur au moins \E9gale \E0 maxChaine
*	lien doit avoir au moins maxMembres elements ( pointeurs sur CROISEMENTS )
*	(11/06/98)	ajout parametre bavard + modif gestion erreurs
*	(16/06/98)	correction pour gerer le noeuds triviaux
*
**********************************************************************/	
int
entrerCodeDynamique(lien,chaine,maxChaine,maxTaille,maxMembres,bavard, p)
CROISEMENT		**lien;
int				maxChaine,maxTaille,maxMembres;
char			*chaine,bavard;
FILE			*p;

{
  char		buff[16],sens,sig,debutMembre=OUI,ttest=NON;
  char *rgs_nb = (char *) NULL,*rgs_si = (char *) NULL,*rgs_pr = (char *) NULL;
  int			l,i,j,k,nbC,maxNum,nbNum,nm,numero;
  CROISEMENT	*debut, *dernier;

  if(chaine==NULL||lien==NULL)return ERROR_RGS_1;
  if(ttest)bavard=OUI;

  /************ allocations dynamiques de memoire ************************/
  maxNum=maxTaille/2;
  if((rgs_nb=malloc(maxNum*sizeof(char)))==NULL){
    if(bavard)report_error (__FILE__, __LINE__,
      "ERREUR : probleme allocation memoire dans entrerCodeDynamique\n");
    return ERROR_RGS_2;
  }
  if((rgs_si=malloc(maxNum*sizeof(char)))==NULL){
    if(bavard)report_error (__FILE__, __LINE__,
      "ERREUR : probleme allocation memoire dans entrerCodeDynamique\n");
    return ERROR_RGS_3;
  }
  if((rgs_pr=malloc(maxNum*sizeof(char)))==NULL){
    if(bavard)report_error (__FILE__, __LINE__,
      "ERREUR : probleme allocation memoire dans entrerCodeDynamique\n");
    return ERROR_RGS_4;
  }
  /************* initialisations tableaux ********************************/
  for(j=0;j<maxNum;j++)rgs_nb[j]=0;
  for(nm=0;nm<maxMembres;nm++)lien[nm]=NULL;

  /************** lecture chaine *****************************************/
  /*if( (fgets(chaine,maxChaine,stdin)==NULL) ){
  return ERREUR_FIN;
  }*/

  //printf("\n%s\n\n", );
  if(fgets(chaine, maxChaine,p) == NULL) return ERREUR_FIN;

  if(ttest)printf("Chaine : (%s)\n",chaine);
  l = strlen(chaine);
  if(l>maxChaine){
    report_error2 (__FILE__, __LINE__,"ATTENTION : Chaine trop grande (max = %d )\n",maxChaine);
    return ERROR_RGS_5;
  }

  /************** traitement *********************************************/
  nm = 0; 											/* numero de membre */
  i = 0;											 /* numero de caractere */
  nbC = 0;										/* numero de croisement */

  /* saute les caracteres de d\E9but *******************/
  while((chaine[i]=='\t' || chaine[i]==' ' 
    || chaine[i]=='e' || chaine[i]=='E')&&(i<l) )i++;

  /**** boucle principale de lecture et de decodage ************************/
  do{
    if(debutMembre){							/* creation debut membre */
      if(ttest)printf("Allocation membre %d\n",nm);
      if((debut = creerCroisement(DEBUT,-1))==NULL){
        printf("ERREUR : allocation memoire\n");
        return ERROR_RGS_6;
      }
      lien[nm]=debut;
      dernier = debut;	
      debutMembre=NON;
    }

    if(ttest)printf("%c(%d)\n",chaine[i],chaine[i]);
    if(chaine[i]=='\n' || chaine[i]=='\r')break;	 /* a cause du gets */

    else if(chaine[i]=='|' ){ 						 /** debut d'un membre **/
      debutMembre=OUI;
      connecteCroisements(dernier,debut);				   /* on boucle */

      nm++;
      if(nm>=maxMembres){
        printf("ERREUR : trop de membres\n");
        return ERROR_RGS_7;
      }
      i++;
    }

    else{
      /*** decode un croisement ******************************************/
      if(chaine[i]=='a') sens = ABOVE;						/* prefixe */
      else  if(chaine[i]=='b') sens = BELOW;
      else {
        if(bavard){
          report_error (__FILE__, __LINE__,"ERREUR : caractere illegal :\n");
          erreurSyntaxe(chaine,i);
        }
        return ERREUR_SYNTAXE;
      }

      if(chaine[i+1]<'0' || chaine[i+1]>'9'){					 /* numero */
        if(bavard){
          report_error (__FILE__, __LINE__,"ERREUR : caractere illegal :\n");
          erreurSyntaxe(chaine,i+1);
        }
        return ERREUR_SYNTAXE;
      }
      j=i+1; k=0;
      while( chaine[j]>='0' && chaine[j]<='9' ) buff[k++]=chaine[j++];
      buff[k]='\0';
      numero=atoi(buff);
      i=j;

#ifdef	CODES_SIGNES
      if(chaine[i]=='+') sig = PLUS;							     /* signe */
      else  if(chaine[i]=='-') sig = MOINS;
      else {
        if(bavard){
          report_error (__FILE__, __LINE__,"ERREUR : caractere illegal :\n");
          erreurSyntaxe(chaine,i);
        }
        return ERREUR_SYNTAXE;
      }
      i++;
#endif

      /**** verifications ************************************************/
      if(numero>maxNum){
        if(bavard){
          report_error2 (__FILE__, __LINE__,"ERREUR : numero trop grand : %d\n",numero);
          erreurSyntaxe(chaine,i+1);
        }
        return ERREUR_SYNTAXE;
      }
      rgs_nb[numero]++;
      if(rgs_nb[numero]>2){
        if(bavard){
          report_error2 (__FILE__, __LINE__,"ERREUR : Trop de passage sur %d\n",numero);
          erreurSyntaxe(chaine,i-2);
        }
        return ERREUR_SYNTAXE;
      }
      if(rgs_nb[numero]==1){						/* 1ere visite   */
        rgs_si[numero]=sig;
        rgs_pr[numero]=sens;
      }
      if(rgs_nb[numero]==2){						/* 2eme visite  */
        if( rgs_si[numero]!=sig ){
          if(bavard){				
            report_error (__FILE__, __LINE__,"ERREUR : Signe incoherent \n");
            erreurSyntaxe(chaine,i-1);
          }
          return ERREUR_SYNTAXE;
        }
        if( rgs_pr[numero]== sens){
          if(bavard){				
            report_error2 (__FILE__, __LINE__,
              "ERREUR : Deux fois meme sens sur %d\n",numero);
            erreurSyntaxe(chaine,i-2);
          }
          return ERREUR_SYNTAXE;
        }
      }


      /*** ajout du croisement dans la liste chainee **********************/ 
      nbC ++;
#ifndef	CODES_SIGNES
      if(ttest)printf("Croisement %d sens %c numero %d\n",nbC+1,
        (sens==ABOVE)?'a':'b', numero); 
      dernier = ajouteCroisement(dernier,numero,sens);
#else
      if(ttest)printf("Croisement %d sens %c numero %d signe %c\n",nbC+1,
        (sens==ABOVE)?'a':'b', numero,(sig==PLUS?'+':'-')); 
      dernier = ajouteCroisement(dernier,numero,sens,sig);
#endif
    }									/** fin boucle decodage croisement **/
  }while(i<l);									 /** fin boucle de lecture **/


  connecteCroisements(dernier,debut);						/* on boucle */
  if(nbC> maxTaille)
    report_error2 (__FILE__, __LINE__,"ATTENTION : Code trop long ( max = %d )\n",maxTaille);

  
  if (rgs_nb) 
    free (rgs_nb);

  if (rgs_si) 
    free (rgs_si);

  if (rgs_pr) 
    free (rgs_pr);

  return nbC;
}


