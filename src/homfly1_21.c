/* fichier homfly1_2l.c
* date creation	:	26/06/97
*
*
* description		:
*	Programme de calcul des polynomes de Homfly pour les liens
* references		:
*	- publication 
*	G. Gouesbet, S. Meunier-Guttin-Cluzel, C. Letellier
*	Computer evaluation of Homfly polynomials by using Gauss codes, with
*	a skein-template algorithm.
*	( soumis a Applied Mathematics and Computations )
*	- rapport technique
*	S. Meunier-Guttin-Cluzel
*	Calcul automatique de polynomes de Homfly associes a des noeuds ou a des
*	liens a partir des codes de Gauss. ( Ref. rtm_07_97_homfly )
*
*
*	Auteur		MEUNIER-GUTTIN-CLUZEL Siegfried
*				LESP-UMR CNRS 6614 Operation 17
*
*	historique		:
*	- Issu de homfly1_1.c
*	- (26/06/97)	gestion des liens
*	- (10/07/97) 	ajout conway et polynomes P
*	version 1.2
*	- 18/06/98	allocations dynamique des tableaux + harmonisation
*
*	fichiers		:
*		croisements_l.c		module de gestion des codes de Gauss
*		croisements_l.o
*		croisements_l.h		
*		polynomesAZ.c		calcul et affichage des polynomes en alpha-z
*		polynomesAZ.h
*		polynomesAZ.o
*		homfly1_2l.c		fichier principal
*
*		
*	cas-test		:	b1-a2-|a3+b4+b2-a1-|a4+b3+           (avec signes)
*		reponse		:   z^-2.(a^-2 -2+a^2) + (a^-2 -2+a^2)-z^2.   w = 0
*
*	compilation cc homfly1_2l.c croisements_l.o polynomesAZ.o -o homl -lm
*		
*	version du 1/07/98
*************************************************************************/

#include		<stdlib.h>
#include		<stdio.h>
#include		<math.h>
#include <string.h>  // RGS
#include		"polynomesAZ.h"

#define			CODES_SIGNES	/** utilise des codes de GAUSS signes */
#include		"cross.h"
#define			BOULE		1
#define			BOULE_PLUS	2
#define			BOULE_MOINS	3
#define			BOULE2		4
#define			SWITCH1		5
#define			SWITCH2		6
#define			SPLICE1		7
#define			SPLICE2		8
#define			A_PLUS		1
#define			A_MOINS		2
#define			B_PLUS		3
#define			B_MOINS		4
#define			HOMFLY_H	0
#define			HOMFLY_P	1
#define			CONWAY		2

#define			TABLE_LISTE		1			/* liste sans numero */
#define			TABLE_NUMERO	2			/* avec un numero */
#define			TABLE_NOEUDS	3			/* 2 numeros */
#define			TABLE_LIENS		4			/* 3 numeros */
#define			TABLE_ETIQUETTE	5			/* une etiquette de type chaine */

#define			ALLOC_INCREMENT	100			/* marge supplementaire pour */
/* augmenter taille polys    */

FILE *fperror = (FILE *) NULL;

void report_error (char *fname, int lnum, char *s) {
  if (!fperror) return;
  fprintf (fperror, "error: line %d in file %s: %s\n", lnum, fname, s);
}
void report_error2 (char *fname, int lnum, char *s, int value) {
  if (!fperror) return;
  fprintf (fperror, "error: line %d in file %s: ", lnum, fname);
  fprintf (fperror, s, value);
  fprintf (fperror, "\n");
}
void report_error3 (char *fname, int lnum, char *s, char *value) {
  if (!fperror) return;
  fprintf (fperror, "error: line %d in file %s: ", lnum, fname);
  fprintf (fperror, s, value);
  fprintf (fperror, "\n");
}



/******************* initialise polynome delta *********************/
void
init_delta(delta)
MONOME			*delta;
{

  delta[0].degA=1;						/** polynome delta **/
  delta[0].degZ=-1;
  delta[0].coeff=1.0;
  delta[1].degA=-1;
  delta[1].degZ=-1;
  delta[1].coeff=-1.0;
}
/******************* affiche message d'erreur memoire ***************/
char
erreur_memoire(bavard,fonction)
char	bavard,*fonction;
{
  if(bavard){
    report_error3 (__FILE__, __LINE__, "ERREUR allocation memoire dans %s\n",fonction);
    return ERREUR_MEMOIRE;
  }
}

short
*rangeNumeros(lien,pNumMin,pNumMax,pNbTriviaux,pNbNum,bavard)
CROISEMENT		**lien;
int				*pNumMin,*pNumMax,*pNbTriviaux,*pNbNum;
char			bavard;
{
  CROISEMENT		*pDeb, *p;
  int				nm,num,tailleNumTab,jMax,j,
    numMin = 32000,numMax = -1,nbTriviaux = 0,NbNum = 0;
  short			*numTab;
  char			ttest=NON;

  if( lien==NULL || lien[0]==NULL ){
    *pNbNum = ERREUR_MEMOIRE;
    return NULL;
  }
  if(ttest)printf("range_numeros entree ------------------\n");
  /*** premier passage : determination min/max et nombre noeuds triviaux ***/
  nm  = 0;
  while( lien[nm]!=NULL ){
    pDeb = saute_debut(lien[nm]);
    if( !est_debut(pDeb) ){
      p = pDeb;
      do{
        num=numero(p);
        if( num>numMax ) numMax = num;
        if( num<numMin ) numMin = num;
        p = suivant_de(p);
      }while( p!=pDeb);
    }
    else{		/** trivial **/
      nbTriviaux++;
    }
    nm++;
  }
  /*** deuxieme passage : allocation de numTab et remplissage *************/
  tailleNumTab = numMax - numMin +1;
  if( ( numTab = (short *)malloc(tailleNumTab*sizeof(short)))==NULL){
    report_error (__FILE__, __LINE__, "ERREUR allocation memoire dans rangeNumeros\n");
    *pNbNum = ERREUR_MEMOIRE;
    return NULL;
  }
  for(j=0;j<tailleNumTab;j++)numTab[j]= -1;

  nm  = 0;
  jMax = 0;
  while( lien[nm]!=NULL ){
    pDeb = saute_debut(lien[nm]);
    if( !est_debut(pDeb) ){
      p = pDeb;
      do{
        num=numero(p);
        j = num - numMin;
        if( numTab[j]== -1 ) numTab[j]=jMax++;
        p = suivant_de(p);
      }while( p!=pDeb);
    }
    nm++;
  }
  *pNumMin = numMin;
  *pNumMax = numMax;
  *pNbTriviaux = nbTriviaux;
  *pNbNum = jMax;
  if(ttest)printf("range_numeros sortie ------------------\n");

  return numTab;
}

/******************* calcul polynomes de HOMFLY, algo STA *************
* version liens (26/06/97)+conway et homfly P (10/06/97)
* (18/06/98)	allocation dynamique + nouvelle gestion des erreurs
**********************************************************************/
int
calculeHomfly(lien,tailleCode,type_polynome,pPolynome,pnPolynome,pWrithe,bavard, pout)
CROISEMENT		**lien;
int				tailleCode;
char			type_polynome,bavard;
MONOME			**pPolynome;
int				*pnPolynome,*pWrithe;
FILE			*pout;
{
  CROISEMENT		*p,*p1,*p2,*pDeb,**pp,**pm;
  int				*nb,*nb1,*deco,
    maxNum,nbNum,numMin,numMax,
    *origine,*final,*numMembre,
    i,ni,i1,j,m,num,
    da,dz,nbMembres,nm,
    nPolt, maxPolt=MAX_POLY, nPol1, maxPol1=MAX_POLY,
    nPolTotal, maxPolTotal=MAX_POLY,
    nPolynome, maxPolynome=MAX_POLY, maxPrevisible,
    nbTriviaux=0, writhe=0;
  short			*numTab;
  char			pref,sig,fini,finBranche,ttest=NON ,branchement,
    queTriviaux=NON;
  MONOME			delta[3],*polt,*pol1,*polTotal,*polynome;
  float			co;

  maxNum=tailleCode/2;

  if( lien==NULL || lien[0]==NULL ) return ERREUR_MEMOIRE;

  if( (numTab = 
    rangeNumeros(lien,&numMin,&numMax,&nbTriviaux,&nbNum,bavard))==NULL )
    return ERREUR_MEMOIRE;
  if(ttest && nbNum!=maxNum)printf("Homfly : ERREUR nbNum incoherent\n");

  /*** allocation de la memoire ***************************************/
  if(ttest)printf("Homfly :alloc. memoire -----------------------------\n");
  if( ( pp= (CROISEMENT **)malloc(tailleCode*sizeof(CROISEMENT)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( ( pm= (CROISEMENT **)malloc(maxNum*sizeof(CROISEMENT)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  /*	if( (numTab = (short *)malloc(nbNum*sizeof(short)))==NULL)
  return erreur_memoire(bavard,"calculeHomfly"); */
  if( (nb= (int *)malloc(nbNum*sizeof(int)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( (nb1= (int *)malloc(tailleCode*sizeof(int)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( (deco= (int *)malloc(nbNum*sizeof(int)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( (origine= (int *)malloc(tailleCode*sizeof(int)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( (final= (int *)malloc(tailleCode*sizeof(int)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( (numMembre= (int *)malloc(tailleCode*sizeof(int)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  /** polynomes ***/
  if( ( polt= (MONOME *)malloc(maxPolt*sizeof(MONOME)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( ( pol1= (MONOME *)malloc(maxPol1*sizeof(MONOME)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( ( polTotal= (MONOME *)malloc(maxPolTotal*sizeof(MONOME)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");
  if( ( polynome= (MONOME *)malloc(maxPolynome*sizeof(MONOME)))==NULL)
    return erreur_memoire(bavard,"calculeHomfly");

  /***************** initialiations **************************************/
  if(ttest)printf("Homfly :init--------------------------------------\n");

  init_delta(delta);

  nPolTotal=0;

  for(j=0;j<maxNum;j++){					 /* initialisations **/
    nb[j]=0; deco[j]=0;
  }
  for(i=0;i<tailleCode;i++)nb1[i]=0;


  i=0;										/* niveau dans l'arbre */
  m=0;						   /* numero membre dans l'etat genere */
  nm=0;						   /* numero membre dans le lien-guide */
  pDeb = saute_triviaux(lien,&nm);		  /* saute noeuds triviaux */
  if(pDeb==NULL && nm==0) return ERREUR_MEMOIRE;
  if(pDeb==NULL) queTriviaux=OUI;
  p=pDeb;
  fini=NON;
  p1=p;										/** debut membre **/
  branchement=NON;


  do{
    finBranche=NON;
    if(!est_debut(p1))do{						/* si pas noeud trivial */
      /*** parcoure membre m *******************************************/
      if(ttest)printf("[%d]",m+1);
      if(!est_debut(p)) do{
        ni=lit_croisement2(p,&pref,&num,&sig,numTab);
        if(ttest)afficheCroisement(stdout,p);

        if(pref==ABOVE)origine[i] = (sig==PLUS ? A_PLUS : A_MOINS );
        else origine[i] = (sig==PLUS ? B_PLUS : B_MOINS );

        pp[i]=p;							 /* enregistre niveau i **/
        numMembre[i]=m;
        nb[num]++;
        nb1[ni]++;


        if( nb[num]==1 ){							/** 1ere visite **/
          if(branchement){
            final[i]=SPLICE1;
            branchement=NON;
          }
          else if( pref==ABOVE ){
            if(sig==PLUS)final[i]=BOULE_PLUS;
            else final[i]=BOULE_MOINS;
          }
          else{
            final[i]= SWITCH1;
          }
          deco[num]=final[i];
        }
        else if( nb[num]==2 ){						/** 2eme visite **/
          if( deco[num]==BOULE_PLUS || deco[num]==BOULE_MOINS)
            final[i]= BOULE2;
          else if( deco[num]==SWITCH1)final[i]= SWITCH2;
          else if( deco[num]==SPLICE1)final[i]= SPLICE2;
        }
        else return erreur_memoire(bavard,"calculeHomfly:visite>2");

        if(ttest){								 /** affichage test **/
          if(final[i]==BOULE2||final[i]==BOULE_PLUS||
            final[i]==BOULE_MOINS)printf("o");
          else if(final[i]==SWITCH1||final[i]==SWITCH2)printf("x");
          else if(final[i]==SPLICE1 ||final[i]==SPLICE2)printf("~");
        }

        if(final[i]==SPLICE1 || final[i]==SPLICE2){		 /* suivant */
          p = suivant_transverse_de(p,lien);
        }
        else p= suivant_de(p);

        i++;
      }while(p!=p1);



      /*** debut membre suivant *****************************************/
      if(ttest)printf("|\n");

      nm = 0;						/* recherche un croisement non visite */
      pDeb = saute_triviaux(lien,&nm);
      if(pDeb!=NULL){
        p2=pDeb;
        ni=lit_croisement2(p2,&pref,&num,&sig,numTab);
        while( nb1[ni]==1 && !finBranche){
          p2 = suivant_de(p2);

          if(p2==pDeb){ 							/* membre suivant */
            nm++;
            pDeb = saute_triviaux(lien,&nm);
            if(pDeb==NULL)finBranche=OUI;
            else p2=pDeb;
          }
          ni=lit_croisement2(p2,&pref,&num,&sig,numTab);
        }
      }
      else finBranche = OUI;
      pm[m]=p1;					/** enregistre debut membre precedent **/

      if(!finBranche){
        p=p2;
        m++;
        p1=p;
      }
    }while(!finBranche);
    if(ttest)printf("\n");

    /*** calculs pour une branche **************/
    da=0; dz=0; co=1.0;
    if(queTriviaux)nbMembres=0; 		/* on a rien fait */
    else nbMembres = m+1;
    nbMembres+= nbTriviaux;
    for(i1=0;i1<tailleCode;i1++){
      if(origine[i1]== A_PLUS && final[i1]==BOULE_PLUS)da++;
      else if(origine[i1]== A_MOINS && final[i1]==BOULE_MOINS)da--;
      else if(origine[i1]==B_PLUS){
        if(final[i1]==SWITCH1)da--;
        else if(final[i1]==SPLICE1)dz++;}
      else if(origine[i1]==B_MOINS){
        if(final[i1]==SWITCH1)da++;
        else if(final[i1]==SPLICE1){
          dz++; co*=-1.0;
        }
      }
    }

    if(ttest)printf("da= %d   dz= %d nm= %d\n",da,dz,nbMembres);
    pol1[0].degA=da;
    pol1[0].degZ=dz;
    pol1[0].coeff=co;
    nPol1=1;
    if(nbMembres>1){
      for(i1=0;i1<nbMembres-1;i1++){
        maxPrevisible = nPol1*2;
        if(maxPrevisible>maxPol1){		/** reallocation des tableaux **/
          maxPol1 = maxPrevisible + ALLOC_INCREMENT;
          if( ( pol1= (MONOME *)
            realloc(pol1,maxPol1*sizeof(MONOME)))==NULL)
            return erreur_memoire(bavard,"calculeHomfly");
          if(bavard)report_error2 (__FILE__, __LINE__, "pol1 passe à %d\n",maxPol1);
          maxPolt = maxPrevisible + ALLOC_INCREMENT;
          if( ( polt= (MONOME *)
            realloc(polt,maxPolt*sizeof(MONOME)))==NULL)
            return erreur_memoire(bavard,"calculeHomfly");
          if(bavard)report_error2 (__FILE__, __LINE__, "polt passe à %d\n",maxPolt);
          if(maxPolt>32000 || maxPol1>32000){
            report_error (__FILE__, __LINE__, "ATTENTION: debordement entier previsible");
            report_error (__FILE__, __LINE__, "dans calculeHomfly\n");
          }
        }
        nPol1=polyProduit(delta,2,pol1,nPol1,polt,nPol1);
        nPol1=polyCopie(pol1,polt,nPol1);
      }
    }
    if(ttest)polyAffiche(pol1,nPol1,pout);

    maxPrevisible = nPol1 + nPolTotal;	/** reallocation si debordement **/
    if(maxPrevisible>maxPolt){
      maxPolt = maxPrevisible;
      if( ( polt= (MONOME *)
        realloc(polt,maxPolt*sizeof(MONOME)))==NULL)
        return erreur_memoire(bavard,"calculeHomfly");
      if(bavard)report_error2 (__FILE__, __LINE__, "polt passe à %d\n",maxPolt);

    }

    nPolt=polySomme(polTotal,nPolTotal,pol1,nPol1,polt,nPolt);

    maxPrevisible = nPolt;				/** reallocation si debordement **/
    if(maxPrevisible>maxPolTotal){
      maxPolTotal = maxPrevisible;
      if( ( polTotal= (MONOME *)
        realloc(polTotal,maxPolTotal*sizeof(MONOME)))==NULL)
        return erreur_memoire(bavard,"calculeHomfly");
      if(bavard)report_error2 (__FILE__, __LINE__, "polTotal passe à %d\n",maxPolTotal);
    }

    nPolTotal=polyCopie(polTotal,polt,nPolt);



    /*** remonte dans l'arbre ********************************************/
    if(i<=0) fini=OUI;
    else{
      i--;
      if(ttest)printf("remonte ");
      do{
        ni=lit_croisement2(pp[i],&pref,&num,&sig,numTab);
        nb[num]--;
        nb1[ni]--;
        i--;
        if(ttest)printf("%d ",i+1);
      }while(i>=0 && final[i]!=SWITCH1);
      if(ttest)printf("\n");
      if(i<0) fini=OUI;
      else{ 									/** SWITCH 1 ***/
        p =pp[i];
        ni=lit_croisement2(p,&pref,&num,&sig,numTab);
        nb[num]--;
        nb1[ni]--;
        m= numMembre[i];
        p1=pm[m];
        branchement=OUI;
      }
    }



  }while(!fini);

  /* calcul du writhe */
  nm=0;
  writhe=0;
  while((pDeb=lien[nm])!=NULL){
    pDeb=saute_debut(pDeb);
    p=pDeb;
    if(!est_debut(pDeb))do{
      writhe += signeCroisement(p);
      p=suivant_de(p);
    }while(p!=pDeb);
    nm++;
  }
  writhe /=2;

  /*** traitement des polynomes speciaux *********************************/
  if(type_polynome == CONWAY){
    for(i=0;i<nPolTotal;i++) polTotal[i].degA=0;
    nPolTotal=polySimplifie(polTotal,nPolTotal);
  }
  if(type_polynome == HOMFLY_P){
    for(i=0;i<nPolTotal;i++) polTotal[i].degA-=writhe;
    nPolTotal=polySimplifie(polTotal,nPolTotal);
  }
  /***********************************************************************/
  if(nPolTotal>maxPolynome){		/** reallocation des tableaux **/
    maxPolynome = nPolynome;
    if( ( polynome= (MONOME *)
      realloc(polynome,maxPolynome*sizeof(MONOME)))==NULL)
      return erreur_memoire(bavard,"calculeHomfly");
    if(bavard)report_error2 (__FILE__, __LINE__, "polynome passe à %d\n",maxPolynome);
  }
  *pnPolynome = polyCopie(polynome,polTotal,nPolTotal);
  *pPolynome = polynome;
  if(ttest)printf("Homfly :retour-------------------------------\n");

  /*** liberation de l'espace memoire ************************************/
  free(pp);
  free(pm);
  free(numTab);
  free(nb);
  free(nb1);
  free(deco);
  free(origine);
  free(final);
  free(numMembre);
  free(polt);
  free(pol1);
  free(polTotal);

  *pWrithe = writhe;
  return 0;
}
/****** Affichage des polynomes : debut de ligne **************************/
void
tableAfficheDebut(option_affichage, type_polynome,tables,n1,n2,n3,etiquette)
char	option_affichage, type_polynome,tables,*etiquette;
int		n1,n2,n3;
{	
  if(option_affichage == LATEX ){
    switch(type_polynome){
      case HOMFLY_H:printf("$H_");break;
      case HOMFLY_P:printf("$P_");break;
      case CONWAY :printf("$C_");break;
    }
    switch(tables){
      case TABLE_NUMERO : 
        //printf("{%d} = $  ",n1);
        break;							
      case TABLE_NOEUDS : 
        //printf("{%d^{%d}} = $  ",n1,n2);
        break;
      case TABLE_LIENS  : 
        //printf("{%d^{%d}_{%d}} = $  ",n1,n2,n3);
        break;
      case TABLE_ETIQUETTE :  
        //printf("{%s} = $  ",etiquette);
        break;
    }
  }
  else if(option_affichage == MAPPLE ){
    switch(type_polynome){
      case HOMFLY_H:printf("H");break;
      case HOMFLY_P:printf("P");break;
      case CONWAY :printf("C");break;
    }
    switch(tables){
      case TABLE_NUMERO : 
        printf("%2.0d := ",n1);break;
      case TABLE_NOEUDS :
        printf("%2.0d%d := ",n1,n2);break;
      case TABLE_LIENS  :
        printf("%2.0d%d%d := ",n1,n2,n3);break;
      case TABLE_ETIQUETTE :
        printf("%s := ",etiquette);break; 
    }
  }
  else{ 
    switch(type_polynome){
      case HOMFLY_H:printf("H");break;
      case HOMFLY_P:printf("P");break;
      case CONWAY :printf("C");break;
    }
    switch(tables){
      case TABLE_NUMERO : 
        printf("(%d) = ",n1);break;							
      case TABLE_NOEUDS : 
        printf("(%d/%d) = ",n1,n2);	break;
      case TABLE_LIENS  : 
        printf("(%d/%d/%d) = ",n1,n2,n3);break;
      case TABLE_ETIQUETTE :  
        printf("(%s) = ",etiquette);break;
    }

  }
}

/*************************  ecran d'aide **********************************/
void
afficheAide()
{
  fprintf(stderr," \n\nHOMFLY1.2         Calcul des polynomes de HOMFLY\n");
  fprintf(stderr," Version liens\n");
  fprintf(stderr," \nMEUNIER-GUTTIN-CLUZEL Siegfried \n\n");
  fprintf(stderr,
    " Syntaxe :   homl [-b][-m][-c xxx][-s xxx][-tx][-ix][-px]\n");
  fprintf(stderr," \n");
  fprintf(stderr," -b    mode bavard\n");
  fprintf(stderr," -m    mode muet\n");
  fprintf(stderr,
    " -c    taille de la chaine de caracteres, defaut -c 1024\n");
  fprintf(stderr," -s    nombre maxi de croisements, defaut -s 256\n");
  fprintf(stderr," -n    nombre maxi de composants, defaut -n 10\n");
  fprintf(stderr," lecture de tables\n");
  fprintf(stderr," -t	table sans numero\n");	
  fprintf(stderr," -tn	table de noeuds (2 numeros)\n");
  fprintf(stderr," -tl	table de liens (3 numeros)\n");
  fprintf(stderr," -t1	table numerotee (1 numero)\n");
  fprintf(stderr," -te	table avec etiquette alphanumerique\n");
  fprintf(stderr," options d'impression (sorties)\n");
  fprintf(stderr," -ie	sortie ecran ( ex : z^2.( 2.a^-2 + 2.a^2) )\n");
  fprintf(stderr," -il	sortie LaTeX \n");
  fprintf(stderr," -im	sortie MAPPLE\n");
  fprintf(stderr," -ib	sortie binaire ( ecrit dans un fichier )\n");
  fprintf(stderr," choix des polynomes \n");
  fprintf(stderr," -ph	polynome H de HOMFLY (defaut)\n");
  fprintf(stderr," -pp	polynome P de HOMFLY \n");
  fprintf(stderr," -pc	polynome de CONWAY\n");
  fprintf(stderr," \n");
  fprintf(stderr," version 18/06/98\n");
}

/**********************************************************************/	
int
main(argc,argv)
int		argc;
char	*argv[];
{
  int kounter = 0;
  MONOME			*polynome;
  CROISEMENT		**lien,	*code;
  int				tailleCode,maxChaine=1024, maxTaille=256,
    maxMembres=MAX_MEMBRES,
    maxPoly=MAX_POLY, nPolynome,maxNum,num,nm,
    i,n1,n2,n3,n4,n5,nn,writhe;
  char			ttest=NON, bavard=NON, detail=NON, 
    valide, tables = NON, filtre=NON,
    cc,*chaine,bid[256],etiquette[32],type_polynome=HOMFLY_P;  // RGS: was HOMFLY_H
  char			input[100], output[100];
  FILE *p = (FILE *) NULL, *pout = (FILE *) NULL;   // RGS
  
  // RGS: input/output file name reading used to be here

  //fperror = fopen ("homfly_errors.txt", "w");
  if (!fperror) fperror = stderr;

  option_affichage = ECRAN;
  i = 1; // RGS
  /*** lecture des parametres de la ligne de commande ********************/
  if(argc>1){
    for(i=1;i<argc;i++){
      if(*(argv[i])=='-'){
        switch((argv[i])[1]){
	case '-':    // RGS
	  p = stdin;
	  pout = stdout;
	  break;
	case 'b' : bavard = OUI; break;
	case 'm' : bavard = NON; break;
	case 'c' :
	  if((i+1)<argc){	
	    maxChaine = atoi(argv[i+1]);
	    if( maxChaine<=0){
	      report_error (__FILE__, __LINE__, "ERREUR : taille illegale\n");
	      afficheAide();
	      exit(-1);
	    }
	    i++;
	  }
	  break;
	case 's' :
	  if((i+1)<argc){	
	    maxTaille = atoi(argv[i+1]);
	    if( maxTaille<=0){
	      report_error (__FILE__, __LINE__, "ERREUR : taille illegale\n");
	      afficheAide();
	      exit(-1);
	    }
	    i++;
	  }
	  break;
	case 'n' :
	  if((i+1)<argc){	
	    maxMembres = atoi(argv[i+1]);
	    if( maxMembres<=0){
	      report_error (__FILE__, __LINE__, "ERREUR : taille illegale\n");
	      afficheAide();
	      exit(-1);
	    }
	    i++;
	  }
	  break;
	case 't' : 				/*** tables de codes *************/
	  switch((argv[i])[2]){
          case '0' :
          case '\0':  tables = TABLE_LISTE; break;
          case '1' :  tables = TABLE_NUMERO; break;
          case '2' :
          case 'n' :  tables = TABLE_NOEUDS; break;
          case '3' :
          case 'l' :  tables = TABLE_LIENS; break;
          case 'e' :  tables = TABLE_ETIQUETTE; break;

          default  : 
            report_error (__FILE__, __LINE__, "ERREUR : option tables illegale\n");
            afficheAide();
            exit(-1);
            }
            break;
	case 'i' :				/*** sortie des polynomes ********/
	  switch((argv[i])[2]){
          case 'l' :  option_affichage = LATEX; break;
          case 'e' :  option_affichage = ECRAN; break;
          case 'm' :  option_affichage = MAPPLE; break;
          case 'b' :  option_affichage = BINAIRE; break;
          default : report_error (__FILE__, __LINE__,
				  "ERREUR : Option affichage illegale\n");
            afficheAide();
            exit(-1);
	  }
	  break;
	  
	case 'p' :				/*** choix des polynomes ********/
	  switch((argv[i])[2]){
          case 'h' : type_polynome=HOMFLY_H; break;
          case 'p' : type_polynome=HOMFLY_P; break;
          case 'c' : type_polynome=CONWAY; break;

          default : report_error (__FILE__, __LINE__,
				  "ERREUR : Choix polynome illegal\n");
            afficheAide();
            exit(-1);
	  }
	  break;
        }
      }
      else break; // RGS
      /*
      else{
        afficheAide();
        exit(-1);
      }
      */
    }
  } 

  // RGS: moved reading input / output files to here

  if (i == argc) {
    if (pout != stdout) {
      printf("input file: "); fflush (stdout);
      scanf("%s", input);
      printf("output file: "); fflush (stdout);
      scanf("%s", output);
    }
  }
  else if (i == argc - 1 && pout != stdout) {
    strcpy (input, argv [i]);
    pout = stdout;
  }
  else if (i == argc - 2 && pout != stdout) {
    strcpy (input, argv [i]);
    strcpy (output, argv [i + 1]);
  }
  else {
    fprintf (stderr, "wrong arguments\n");
    exit (102);
  }

  if (p != stdin) {
    p=fopen(input, "r");
    if (!p) {
      fprintf (stderr, "Can't open file `%s' for reading\n", input);
      exit (102);
    }
  }
  if (pout != stdout) {
    pout=fopen(output, "w+");
    if (!pout) {
      fprintf (stderr, "Can't open file `%s' for output\n", output);
      exit (102);
    }
  }
  /*************************** Allocation memoire ***********************/

  if((chaine=malloc(maxChaine*sizeof(char)))==NULL){
    report_error (__FILE__, __LINE__, "ERREUR : probleme allocation memoire dans main\n");
    exit(-1);
  }
  if((polynome=(MONOME *)malloc(maxPoly*sizeof(MONOME)))==NULL){
    report_error (__FILE__, __LINE__, "ERREUR : probleme allocation memoire dans main\n");
    exit(-1);
  }
  if((lien=(CROISEMENT **)malloc(maxMembres*sizeof(CROISEMENT *)))==NULL){
    report_error (__FILE__, __LINE__, "ERREUR : probleme allocation memoire dans main\n");
    exit(-1);
  }

  /*** entree du code **************************************************/
  //fprintf(stderr,"Calcul automatique des polynomes de HOMFLY\n");
  //fprintf(stderr,"	a partir des codes de GAUSS signes pour des liens.\n");
  //fprintf(stderr,"Entrer le code signe du lien (a1+b2-|a3+...)\n");
  //fprintf(stderr," ou controle-D pour arreter\n");



  /*** boucle de lecture et de traitement *******************************/
  //int index;

  do{
    nn = 1;
    if(tables>0){
      /** saute les commentaires ***/
      while((cc=getc(stdin))=='#'||cc=='%')gets(bid);
      ungetc(cc,stdin);
      /** lit les numeros ***/
      switch(tables){
        case TABLE_NUMERO : nn=scanf("%d ",&n1); break;
        case TABLE_NOEUDS : nn=scanf("%d %d ",&n1,&n2); break;
        case TABLE_LIENS : nn=scanf("%d %d %d ",&n1,&n2,&n3); break;
        case TABLE_ETIQUETTE : nn=scanf("%s ",etiquette); break;
      }
    }
    //else fprintf(stderr,"code ? :");

    if(nn<1) tailleCode=ERREUR_FIN;
    else {
      /*** entree et test *******************************************/

      tailleCode = entrerCodeDynamique
        (lien,chaine,maxChaine,maxTaille,maxMembres,bavard,p);

      //fprintf (pout, "%d ", kounter);
      ++kounter;
      if(tailleCode>0){ 

          if(ttest)afficheCode(lien);
          /*** affichage debut de ligne *****************************/
          if(tables>0)tableAfficheDebut
            (option_affichage,type_polynome,tables,n1,n2,n3,etiquette);

          /************** calcul et sortie polynome *****************/
          if(calculeHomfly(lien,tailleCode,type_polynome,
            &polynome,&nPolynome,&writhe,bavard) == ERREUR_MEMOIRE) {
              fprintf (fperror, "ERREUR_MEMOIRE occured\n");
              exit(-1);
          }

          if(option_affichage == LATEX)
            polyLatex(polynome,nPolynome);
          else if(option_affichage == MAPPLE)
            polyMapple(polynome,nPolynome);
          else
            polyAffiche(polynome,nPolynome, pout);


          /*** affichage fin de ligne *****************************/
          if(option_affichage == LATEX )
            printf(" $  [\\omega = %d]$ \\\\\n",writhe);
          /*else if(option_affichage != MAPPLE)
          printf(" w = %d\n\n",writhe);*/
      }
      else {
        if (fperror && tailleCode != -3) 
          fprintf (fperror, "on line %d error code was %d\n", kounter, tailleCode);
      }

      nm=0;							/** liberation de la memoire **/
      while(lien[0]!=NULL){
        libereCode(lien,0);
      }
      if(ttest)printf("Taille Code : %d\n",tailleCode);
    }
  }while(tailleCode>1);

  if (p != stdin) fclose(p);       // RGS
  if (p != stdout) fclose(pout);   // RGS
  if (fperror && fperror != stderr) fclose (fperror);

  return 0;	
}
