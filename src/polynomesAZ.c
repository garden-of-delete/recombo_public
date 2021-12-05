/* fichier	polynomesAZ.c
 * date creation	:	20/06/97
 *	
 * description		:
 *
 *	Fonctions de manipulation de polynomes en a-z, utilises pour les
 *	calculs d'invariants polynomiaux pour les liens et les noeuds
 *
 *	Auteur		MEUNIER-GUTTIN-CLUZEL Siegfried
 *				LESP-UMR 6614 Coria- Operation 17
 *
 *	historique		:
 *						- 22/10/97		amelioration pour degre alpha =0
 *						- 27/11/97		correction pour alpha=0 dz=0
 *						- 14/05/98		ajout sortie compatible MAPPLE
 *										allocation dynamique
 *
 *	compilation		:	cc -c polynomesAZ.c
 *	version du 14/05/98
 ***********************************************************************/
#include		<stdlib.h>
#include		<stdio.h>
/*#include		<math.h>*/
#include		"polynomesAZ.h"


/********************* allocation dynamique ****************************/
int
polyAlloue(pp1,n1)
	MONOME		**pp1;
	int			n1;
{
	if(pp1==NULL || n1<1)return -1;

	*pp1 = (MONOME *)malloc(n1*sizeof(MONOME));
	if(*pp1==NULL) return -1;
	else return n1;

}
int
polyDesalloue(p1)
	MONOME		*p1;
{
	if(p1==NULL)return -1;

	free(p1);
	return 0;
}

/********************* copie p2 dans p1 retourne n2 ********************/
int
polyCopie(p1,p2,n2)
	MONOME		p1[],p2[];
	int			n2;
{
	int			i1;

	for(i1=0;i1<n2;i1++){
		p1[i1].degA = p2[i1].degA;
		p1[i1].degZ = p2[i1].degZ;
		p1[i1].coeff =p2[i1].coeff;
	}
	return n2;
}


int
polyProduit(p1,n1,p2,n2,p3,n3)
	MONOME		p1[],p2[],p3[];
	int			n1,n2,n3;
{
	int			i1,i2;

/*	printf("polyproduit : %d, %d\n",n1,n2); */
	for(i1=0;i1<n1;i1++)
		for(i2=0;i2<n2;i2++){
			p3[i1*n2+i2].coeff= (p1[i1].coeff)*(p2[i2].coeff);
			p3[i1*n2+i2].degA = p1[i1].degA+ p2[i2].degA;
			p3[i1*n2+i2].degZ= p1[i1].degZ + p2[i2].degZ;
		}
	n3=n1*n2;
	return polySimplifie(p3,n3);
}

void
polyMultiplie(p1,n1,da,dz,c)
	MONOME		p1[];
	int			n1,da,dz;
	float		c;
{
	int			i1;

	for(i1=0;i1<n1;i1++){
		p1[i1].degA+= da;
		p1[i1].degZ+= dz;
		p1[i1].coeff*=c;
	}
}
	
int
polySomme(p1,n1,p2,n2,p3,n3)
	MONOME		p1[],p2[],p3[];
	int			n1,n2,n3;
{
	int			i1;

/*	printf("polysomme : %d, %d\n",n1,n2); */
	for(i1=0;i1<n1;i1++){
		p3[i1].coeff = p1[i1].coeff ;
		p3[i1].degA = p1[i1].degA ;
		p3[i1].degZ = p1[i1].degZ ;
	}
	for(i1=0;i1<n2;i1++){
		p3[i1+n1].coeff = p2[i1].coeff ;
		p3[i1+n1].degA = p2[i1].degA ;
		p3[i1+n1].degZ = p2[i1].degZ ;
	}
	n3=n1+n2;
	return polySimplifie(p3,n3);
}

int
polySimplifie(p1,n1)
	MONOME		p1[];
	int			n1;
{
	int			i1,i2;

	for(i1=0;i1<n1-1;i1++)
		for(i2=i1+1;i2<n1;i2++){
			if((p1[i1].degA==p1[i2].degA)&&( p1[i1].degZ==p1[i2].degZ)){
				p1[i1].coeff+=p1[i2].coeff;
				p1[i2].coeff=0;
			}
	}
	/** elimination des elements nuls **/
	i1=0; i2=0;
	while(i1<n1){

		if(i1>i2 && p1[i1].coeff!=0){
			p1[i2].coeff=p1[i1].coeff;		
			p1[i2].degZ=p1[i1].degZ;
			p1[i2].degA=p1[i1].degA;
			i2++;
		}
		else if( i1==i2 && p1[i1].coeff!=0)i2++;
		i1++;
	}
	return i2;
}

void polyTrie(p1,n1)
	MONOME		p1[];
	int			n1;
{
	int			i1,i2,it;
	float		ft;

	for(i1=0;i1<n1-1;i1++){
		for(i2=i1+1;i2<n1;i2++){
			if( (p1[i1].degZ>p1[i2].degZ)||
				((p1[i1].degZ==p1[i2].degZ)&&(p1[i1].degA>p1[i2].degA))){
				it=p1[i1].degZ; p1[i1].degZ=p1[i2].degZ; p1[i2].degZ=it;
				it=p1[i1].degA; p1[i1].degA=p1[i2].degA; p1[i2].degA=it;
				ft=p1[i1].coeff; p1[i1].coeff=p1[i2].coeff; p1[i2].coeff=ft;
			}
		}
	}
}
	
/***** Affichage du polynome a l'ecran ************************************/
void polyAffiche(p1,n1, pout)
	MONOME		p1[];
	int			n1;
	FILE		*pout;
{
	int			i1,nz,i,dz,da;
	float		co;

	polyTrie(p1,n1);				/* trie les monomes */


	if(n1==0){						/* polynome nul */
		fprintf(pout, "0\n");
		return;
	}
	i1=0;
	while(i1<n1){
		dz=p1[i1].degZ;					/* compte les monomes de degre dz */
		nz=1;
		i=i1+1;
		while(p1[i].degZ==dz && i++<n1)nz++;

		if(nz==1){							/* pas besoin de parentheses */
			co = p1[i1].coeff;
			da = p1[i1].degA;
			/*** Co ***********************/
			if(i1>0 && co>0.0)fprintf(pout, " + ");
			if(co!=1 && co!=-1){
				fprintf(pout, "%.0f",co);
				if(dz!=0 || da!=0) fprintf(pout, ".");
			}
			else if(co==-1)fprintf(pout, " -");
			if(dz==0 && da==0)fprintf(pout, "1");
			/*** z ***********************/
			if(dz==1)fprintf(pout, "z");
			else if(dz!=0)fprintf(pout, "z^%d",dz);
			if(dz!=0 && da!=0) fprintf(pout, ".");
			/*** Alpha ***********************/
			if(da==1)fprintf(pout, "a");
			else if(da!=0)fprintf(pout, "a^%d",da);



		}
		else{									/* on a une parenthese */	
			if(i1>0)fprintf(pout, " + ");
		
			/*** z ***********************/
			if(dz==1)fprintf(pout, "z.");
			else if(dz!=0)fprintf(pout, "z^%d.",dz);
		
			fprintf(pout, "( ");
			for(i=i1;i<i1+nz;i++){			  /* boucle sur les degres de z */
				co = p1[i].coeff;
				da=p1[i].degA;
				/*** Co ***********************/
				if(i>i1 && co>0.0)fprintf(pout, " +");
				if(da==0){
					fprintf(pout, " %.0f",co);
				}
				else{
					if(co!=1 && co!=-1){
						fprintf(pout, " %.0f",co);
						if(da!=0)fprintf(pout, ".");
					}
					else if(co==-1)fprintf(pout, "-");
				}	
				/*** Alpha ***********************/
				if(da==1)fprintf(pout, "a");
				else if(da!=0)fprintf(pout, "a^%d",da);
			}
			fprintf(pout, " )");
		}
		i1+=nz;
	}
	fprintf(pout, "\n");
}
/***** Affichage du polynome compatible LaTeX *****************************/
void polyLatex(p1,n1)
	MONOME		p1[];
	int			n1;
{
	int			i1,nz,i,dz,da;
	float		co;

	polyTrie(p1,n1);								/* trie les monomes */

	printf("$ ");
	if(n1==0){											/* polynome nul */
		printf("0$\n");
		return;
	}
	i1=0;
	while(i1<n1){
		dz=p1[i1].degZ;			  /* compte les monomes de degre dz en z */
		nz=1;
		i=i1+1;
		while(p1[i].degZ==dz && i++<n1)nz++;

		if(nz==1){							/* pas besoin de parentheses */
			co = p1[i1].coeff;
			da = p1[i1].degA;
			/*** Co ***********************/
			if(i1>0 && co>0.0)printf(" + ");
			if(co!=1 && co!=-1){
				printf("%.0f",co);
				if(dz!=0 || da!=0) printf("*");
			}
			else if(co==-1)printf(" -");
			if(dz==0 && da==0)printf("1");
			/*** z ***********************/
			if(dz==1)printf("z");
			else if(dz!=0)printf("z^{%d}",dz);
			if(dz!=0 && da!=0) printf(".");
			/*** Alpha ***********************/
			if(da==1)printf("\\alpha");
			else if(da!=0)printf("\\alpha^{%d}",da);

		}
		else{									/* on a une parenthese */		
			if(i1>0)printf(" + ");
		
			/*** z ***********************/
			if(dz==1)printf("z.");
			else if(dz!=0)printf("z^{%d}.",dz);
		
			printf("( ");
			for(i=i1;i<i1+nz;i++){			  /* boucle sur les degres de z */
				co = p1[i].coeff;
				da = p1[i].degA;
				/*** Co ***********************/
				if(i>i1 && co>0.0)printf(" +");
				if(da==0){
					printf(" %.0f",co);
				}
				else{
					if(co!=1 && co!=-1){
						printf(" %.0f",co);
						if(da!=0)printf(".");
					}
					else if(co==-1)printf("-");
				}	
				/*** Alpha ***********************/
				if(da==1)printf("\\alpha");
				else if(da!=0)printf("\\alpha^{%d}",da);
			}
			printf(" )");
		}
		i1+=nz;
	}
	printf("$\n");
}
/***** Affichage du polynome compatible MAPPLE ********(14/05/98)***********/
void polyMapple(p1,n1)
	MONOME		p1[];
	int			n1;
{
	int			i1,nz,i,dz,da;
	float		co;

	polyTrie(p1,n1);				/* trie les monomes */


	if(n1==0){						/* polynome nul */
		printf("0;\n");
		return;
	}
	i1=0;
	while(i1<n1){
		dz=p1[i1].degZ;				/* compte les monomes de degre dz en z */
		nz=1;
		i=i1+1;
		while(p1[i].degZ==dz && i++<n1)nz++;

		if(nz==1){					/* un seul monome => pas de parentheses */
			co = p1[i1].coeff;
			da = p1[i1].degA;
			/*** Co ***********************/
			if(i1>0 && co>0.0)printf(" + ");
			if(co!=1 && co!=-1){
				printf("%.0f",co);
				if(dz!=0 || da!=0) printf("*");
			}
			else if(co==-1)printf(" -");
			if(dz==0 && da==0)printf("1");
			/*** z ***********************/
			if(dz==1)printf("z");
			else if(dz!=0)printf("z**%d",dz);
			if(dz!=0 && da!=0) printf("*");
			/*** Alpha ***********************/
			if(da==1)printf("a");
			else if(da!=0)printf("a**%d",da);
		}
		else{							/* plusieurs monomes => parentheses */					
			if(i1>0)printf(" + ");
		
			/*** z ***********************/
			if(dz==1)printf("z*");
			else if(dz!=0)printf("z**%d*",dz);
		
			printf("( ");
			for(i=i1;i<i1+nz;i++){			  /* boucle sur les degres de z */
				co = p1[i].coeff;
				da = p1[i].degA;
				/*** Co ***********************/
				if(i>i1 && co>0.0)printf(" +");
				if(da==0){
					printf(" %.0f",co);
				}
				else{
					if(co!=1 && co!=-1){
						printf(" %.0f",co);
						if(da!=0)printf("*");
					}
					else if(co==-1)printf("-");
				}	
				/*** Alpha ***********************/
				if(da==1)printf("a");
				else if(da!=0)printf("a**%d",da);
			}
			printf(" )");
		}
		i1+=nz;
	}
	printf(";\n");
}

/************* entree manuelle d'un polynome *******************************/

int
polyEntre(p1,n1)
	MONOME		p1[];
	int			n1;
{
	int			i1;

	i1=0;
	do{
		printf("Coeff  : ");
		scanf("%f",&(p1[i1].coeff));
		if(p1[i1].coeff!=0){
			printf("deg A  : ");
			scanf("%d",&(p1[i1].degA));
			printf("deg Z  : ");
			scanf("%d",&(p1[i1].degZ));
		}
		i1++;
	}
	while(i1<MAX_POLY &&  p1[i1-1].coeff!=0);

	return i1-1;
}
