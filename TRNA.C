/*****************************************************************************
* RNA structure secondaire par glissement                                    *
*                                                                            *
* Ce programme est la version finale du programme RNA.  Il reconnait les     *
* bases suivantes:                                                           *
*       A = Adenine             C = Cytidine                                 *
*       G = Guanosine           B = Bromouridine                             *
*       I = Inosine             D = Dihydrouridine                           *
*       M = Thioinosine         P = Pseudouridine                            *
*       X = Xanthosine          T = Ribothymidine                            *
*       R = Une Puranosine      S = Thiouridine                              *
*       N = Un nucleoside       U = Uridine                                  *
*                               Y = Une Pyrimididine                         *
*                                                                            *
* De plus, il accepte un seul type de ponts avec trois liens, et ils doivent *
* etre composes seulement de liens C-G. Les boucles plus petite que la       *
* taille minimale entree par l'utilisateur sont eliminees.                   *
*                                                                            *
* Les kgraphes contenant au moins series de 4 ou 5 bons ponts (branches en   *
* helice) sont seulement retenus. Ceux de 6 ou plus et ceux de 3 ou moins    *
* sont rejetes.                                                              *
*                                                                            *
*    Les fichiers crees par ce programme sont les suivants:                  *
*                                                                            *
*    fg... : contient chaque etape du glissement de la molecule.             *
*    fq... : contient les bons ponts du fichier fg...                        *
*    fm... : matrice des ponts.                                              *
*    fc... : matrice de compatibilite.                                       *
*    fk... : l'identification des kgraphes.                                  *
*    fl... : l'identification avec les lignes correspondantes                *
*            de la matrice des ponts pour chaque kgraphe.                    *
*    fd... : l'identification, les liens, le nombre de ponts,                *
*            la grosseur et le type des boucles et l'energie des ponts.      *
*    fe... : l'identification et le dessin de chaque kgraphe.                *
*    ft... : tableau final.                                                  *
*    fx... : le numero des lignes de la matrice des ponts avec               *
*             les liens correspondants.                                      *
*    ff... : diagrammes de connectivite                                      *
*    fp... : compte des kgraphes et des dessins                              *
*                                                                            *
*****************************************************************************/


#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <string.h>

#define NEWTABLE(var,quan) var = (void *) malloc(quan*sizeof(var[0]))
#define AND &&
#define OR ||

FILE *in;
FILE *in1;
FILE *out;
FILE *out1;
FILE *out2;

int longue, dim=0;
int knombre=0;
int k4count=0,k5count=0;
static int counts[5][19] = {
				 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
			   };
int nboucle=0;
int ncouple=0;
int petite;
int pl[21];
static int pp[7];
int flag2=0;
char *chaine, *ch_inv;
char name[4];
char rep,rep2;
static char gu[8]="G";
static char cr[8]="C";
static char fg[8]="FG";
static char fq[8]="FQ";
static char fc[8]="FC";
static char fm[8]="FM";
static char fl[8]="FL";
static char fk[8]="FK";
static char fd[8]="FD";
static char ft[8]="FT";
static char fe[8]="FE";
static char fx[8]="FX";
static char ff[8]="FF";
static char fp[8]="FP";
		
void fich_fg();
void fich_fgs();
void pont1();
void pont2();
void mat_loop();
int boucle_trop_petite();
void mat1();
void mat2();
void croise();
void kgraphe();
void marque_pont();
void plus_petit();
void couplage();
void calcul_pont_boucle();
void dessin();
void ecrire_fd();
int retrouve_pont_interieur();
float calcul_energie();
float nombre_boucle();


void trna1()
/*****************************************************************************
* Lecture du fichier RNA contenant la chaine a etre etudiee                  *
*****************************************************************************/
  {
	static char rna[7]="RNA";
	int temp;
	int indice = 0;

	printf("\n\n   *****   tRNA Secondary structure prediction program   *****\n\n");
	do
	  {
	printf(" tRNA primary structure (sequence) file: RNA");
	scanf("%3s",name);
	  }
	while(strlen(name)==0);
	strcat(rna,name);
	do
	  {
	printf("Do you wish to consider non-Watson-Crick cases [Y / N] ? ");
	rep=(char)getch();
	printf("%c\n",rep);
	  }
	while(rep!='Y' AND rep!='y' AND rep!='N' AND rep!='n');
	do
	  {
	printf("Do you wish to consider knotted structures [Y / N] ? ");
	rep2=(char)getch();
	printf("%c\n",rep2);
	  }
	while(rep2!='Y' AND rep2!='y' AND rep2!='N' AND rep2!='n');
	do
	  {
	printf("What is the smallest loop size to be considered ? ");
	scanf("%d",petite);
	  }
	while(petite<0);
	if ( (in = fopen(rna,"r") ) != NULL)
	  {
	indice=0;
	while(!feof(in)) for( ; isalnum(temp=getc(in)); indice++);
	longue=indice;
	rewind(in);                /* Retourne pointeur au debut du fichier */
	chaine = (char *) malloc(longue +2);
	if (chaine == NULL)
	  {
		printf("Allocation error - end");
		exit(2);
	  }
	indice=0;
	while(!feof(in)) for( ; isalnum(temp=getc(in)); chaine[indice++]=(char)temp);
	chaine[indice] = '\0';
	fclose(in);
	  }
	  else
	{
	  printf("The tRNA sequence input file %s cannot be opened.\n",rna);
	  exit(1);
	}
  }
/*****************************************************************************
* Fin de la fonction trna1()                                                 *
*****************************************************************************/


void init_fichiers()
/*****************************************************************************
* fonction qui vide les fichiers a etre ecrits                               *
*****************************************************************************/
  {
	static char name2[5];

	if((rep=='y' OR rep=='Y') AND (rep2=='n' OR rep2=='N'))
	  {
	strcat(gu,name);
	strcpy(name2,gu);
	  }
	  else if((rep=='n' OR rep=='N') AND (rep2=='y' OR rep2=='Y'))
	{
	  strcat(cr,name);
	  strcpy(name2,cr);
	}
	else if((rep=='y' OR rep=='Y') AND (rep2=='y' OR rep2=='Y'))
	  {
		strcat(gu,cr);
		strcat(gu,name);
		strcpy(name2,gu);
	  }
	  else
		strcpy(name2,name);
	strcat(fg,name2);
	strcat(fq,name2);
	strcat(fm,name2);
	strcat(fc,name2);
	strcat(fk,name2);
	strcat(fl,name2);
	strcat(fd,name2);
	strcat(ft,name2);
	strcat(fe,name2);
	strcat(fx,name2);
	strcat(ff,name2);
	strcat(fp,name2);
	if ((out = fopen(fg,"w")) != NULL) fclose(out);
	if ((out = fopen(fq,"w")) != NULL) fclose(out);
	if ((out = fopen(fm,"w")) != NULL) fclose(out);
	if ((out = fopen(fc,"w")) != NULL) fclose(out);
	if ((out = fopen(fk,"w")) != NULL) fclose(out);
	if ((out = fopen(fl,"w")) != NULL) fclose(out);
	if ((out = fopen(fd,"w")) != NULL) fclose(out);
	if ((out = fopen(ft,"w")) != NULL) fclose(out);
	if ((out = fopen(fe,"w")) != NULL) fclose(out);
	if ((out = fopen(fx,"w")) != NULL) fclose(out);
	if ((out = fopen(ff,"w")) != NULL) fclose(out);
	if ((out = fopen(fp,"w")) != NULL) fclose(out);
  }
/*****************************************************************************
* Fin de la fonction init_fichiers()                                         *
*****************************************************************************/


void invers()
/*****************************************************************************
* Fonction d'inversion de la chaine de RNA                                   *
*****************************************************************************/
  {
	int indice;
	char *indexe, *ind_inv, *addre1, *addre2;

	ch_inv = (char *) malloc(longue +2);
	if (ch_inv==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	indexe = (char *) malloc(longue +2);
	if (indexe==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	ind_inv = (char *) malloc(longue +2);
	if (ind_inv==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	for ( indice = 1; indice <= longue ; indice++)
	  {
	indexe[indice] = (char)(indice % 10);
	printf("%d",indexe[indice]);
	  }
	printf("\n");
	addre1 = chaine + longue;
	addre2 = ch_inv;
	while (--addre1 >= chaine) *addre2++ = *addre1;
	ch_inv[longue] = '\0';
	addre1 = indexe +1 + longue;
	addre2 = ind_inv;
	while (--addre1 >= indexe) *addre2++ = *addre1;
	ind_inv[longue] = '\0';
	printf("%s\n%s\n",chaine, ch_inv);
	for(indice=0;indice<longue; indice++)
	  printf("%d",ind_inv[indice]);
	printf("\n");
	if ((in = fopen(fg,"a")) != NULL)
	  {
	fich_fg(indexe,longue+1,1,1);
	fich_fg(chaine,longue,0,0);
	fprintf(in,"\n");
	fich_fg(ch_inv,longue,0,0);
	fich_fg(ind_inv,longue,0,1);
	fprintf(in,"\n\n");
	fich_fg(indexe,longue+1,1,1);
	fclose(in);
	  }
	  else
	{
	  printf("The file %s could not be opened.\n",fg);
	  exit(3);
	}
	free(indexe);
	free(ind_inv);
/*****************************************************************************
* Remet la memoire utilisee par indexe et ind_inv disponible.                *
*****************************************************************************/
  }
/*****************************************************************************
* Fin de la fonction invers()                                                *
*****************************************************************************/


void fich_fg(char *chaine1,int long3, int indice, int jj)
/*****************************************************************************
* Fonction permettant d'ecrire des chaines de caracteres ou d'entiers.       *
* Definition des parametres:                                                 *
*       chaine1: chaine a etre tranferee dans le fichier fg.                 *
*       long3:   borne superieure qui n'est pas inclut dans l'impression.    *
*       indice:  borne inferieure d'indexation de la chaine.                 *
*       jj:      si != 1 > imprimera des caracteres, sinon imprimera des     *
*                entiers.                                                    *
*****************************************************************************/

  {
	static int ii;

	if(jj!=1)
	  {
	for (ii = indice;ii < long3; ii++) fprintf(in,"%c",chaine1[ii]);
	fprintf(in,"\n");
	  }
	  else
/*****************************************************************************
* Si j=1: alors chaine1 contient des entiers                                 *
*****************************************************************************/
	{
	  for (ii = indice;ii < long3; ii++) fprintf(in,"%d",chaine1[ii]);
	  fprintf(in,"\n");
	}
  }
/*****************************************************************************
*  Fin de la fonction fich_fg()                                              *
*****************************************************************************/


void glisse()
/*****************************************************************************
* Fonction de glissement de la molecule sur son inverse.  Et definition des  *
* symboles.                                                                  *
*****************************************************************************/
  {
	static char *ch_bond;
	static int ind,i,j,g,compteur,longx2,l,k,k1;
	static char bp[3];
	static char possA[] = "AG.GA.AI.IA.AM.MA.AX.XA.BC.CB.BG.GB.BI.IB.BM.MB.CD.DC.CO.OC.CS.SC.CT.TC.CU.UC.CP.PC.DG.GD.DI.ID.DM.MD.GO.OG.GS.SG.GT.TG.GU.UG.GP.PG.IO.OI.IU.UI.IS.SI.IT.TI.IP.PI.MO.OM.MS.SM.MT.TM.MU.UM.MP.PM" ;
	static char possB[] = "AB.BA.AD.DA.AO.OA.AS.SA.AT.TA.AU.UA.AP.PA.GI.IG.GM.MG.GX.XG" ;
	static char possC[] = "CG.GC" ;

	ch_bond = (char *) malloc(longue +2);
	if (ch_bond==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	strcpy(ch_bond,"");
	longx2 = 2*longue;
	printf("\n\n");
	compteur=0;
	ind = 1;
	if ((in = fopen(fg,"a")) != NULL)
	  {  
	bp[2] = '\0';
	for( i=longx2; i>1; i--)
	  {
		k1 = (i>longue) ? i - longue-1:0;
		g = (i <= longue) ?  ind++ :0;
		if (k1 != 0) l = (longue - k1)/2;
		  else l = (longue + g)/2;
		for (j = g, k = k1;j < l AND k < longue; j++, k++)
		  {
		bp[0] = chaine[k];
		bp[1] = ch_inv[j];
		if(strstr(possA,bp) != NULL )  strcat(ch_bond,"-");
		  else if(strstr(possB,bp) != NULL ) strcat(ch_bond,"+");
			else if(strstr(possC,bp) != NULL ) strcat(ch_bond,"*");
			  else strcat(ch_bond,".");
		  }
/*****************************************************************************
* Fin de la boucle j                                                         *
*****************************************************************************/
		++compteur;
		fich_fgs(compteur,ch_bond,k1+1, longue -g);
		strcpy(ch_bond,"");
	  }
/*****************************************************************************
* fin de la boucle for i                                                     *
*****************************************************************************/
	fclose(in);
	  }
	  else
	{
	  printf("The file %s could not be opened.\n",fg);
	  exit(3);
	}
	free(ch_bond);
	free(ch_inv);
/*****************************************************************************
* Libere la memoire utilisee par ces variables disponible pour la fonction      *
* malloc()                                                                   *
*****************************************************************************/
  }
/*****************************************************************************
* Fin de la fonction glisse()                                                *
*****************************************************************************/


void fich_fgs(int compt,char *string,int ii,int jj)
/*****************************************************************************
* ajoute les strings de bonds au fichier fg...                               *
*****************************************************************************/
  {
	if(strlen(string)>2) 
	  {
	fprintf(in,"^%d>     %d - %d\n",compt,ii,jj);
	fprintf(in,"%s\n",string);
	  }
  }
/*****************************************************************************
*   fin de la fonction fich_fgs()                                            *
*****************************************************************************/


void cherche_pont()
/*****************************************************************************
* Fonction de recherche des ponts: des differentes combinaisons possibles de *
* * et + et des cas ou le - n'est pas a l'une des extremites                 *
*****************************************************************************/
  {
	char *pont;
	int i;

	pont = (char * ) malloc(longue+2);
	if (pont==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	for( i=0;i<longue+1;i++) pont[i]='\0';
	if((in=fopen(fg,"r"))==NULL)
	  {
	printf("Error opening file %s for reading - end\n",fg);
	exit(1);
	  }
	if((out=fopen(fq,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",fq);
	exit(3);
	  }
	if (rep == 'n' OR rep=='N') pont1(pont);
	   else pont2(pont);
	fclose(in);
	fclose(out);
	free(pont);
  }
/*****************************************************************************
* fin de la fonction cherche_pont()                                          *
*****************************************************************************/


void pont1(char *pont)
/*****************************************************************************
* sans l'option G - U.                                                       *
* Recherche des ponts ****, **-**, ***-** (ou les etoiles peuvent etre       *
* remplacees par des plus).  C.-a-d. sans l'option G-U                       *
*****************************************************************************/
  {
	static char temp;
	int i,j,i1,i2,i3,len;
	static int flag3=0;

	while( (int)(temp=(char)getc(in)) !=EOF)
	  {
	while((int)temp != EOF AND (int)temp != 94) temp = (char)getc(in);
	fscanf(in,"%d>     %d - %d\n",&i1,&i2,&i3);
	len=(i3-i2)/2+(i3-i2)%2;
	if(fgets(pont,len+1,in) != NULL AND i1>5 AND i3>5)
	  {
		pont[len+1]='\0';
		len=strlen(pont);
		for(i=0;i<len-2;i++)
		  {
		if((pont[i]=='+' OR pont[i]=='*')AND(pont[i+1]=='*'OR pont[i+1]=='+')AND pont[i+2] != '.' )
		  {
			if(pont[i+2]=='-' AND ((pont[i+3]!='*'AND pont[i+3]!='+') OR (pont[i+4]!='+' AND pont[i+4]!='*')));
			  else
			if((pont[i]!='*' OR pont[i+1]!='*' OR pont[i+2]!='*') AND (pont[i+3]!='-' AND pont[i+3]!='+' AND pont[i+3]!='*'));
			  else
				if((pont[i]!='*' OR pont[i+1]!='*' OR pont[i+2]!='*') AND pont[i+3]=='-' AND ((pont[i+4]!='+' AND pont[i+4]!='*') OR (pont[i+5]!='+' AND pont[i+5]!='*')));
				  else
				{
				  j = i+4;
				  while(pont[j]=='+' OR pont[j]=='*')
					{
					  if(i2+j == i3-j+1)
					{
					  flag3=1;
					  break;
					}
					  j++;
					}
				  if(flag3==0)
					{
					  fprintf(out,"%d          %d   %d - %d   %d\n",i1,i,i2,i3,len);
					  fprintf(out,"%s\n",pont);
					}
				}
			break;
		  }
		  }
/*****************************************************************************
* fin de la boucle for i                                                     *
*****************************************************************************/
	  }
	  }
  }
/*****************************************************************************
* fin de l'ensemble sans l'option G-U                                        *
*****************************************************************************/


void pont2(char *pont)
/*****************************************************************************
* debut de l'ensemble avec l'option G - U.                                   *
*****************************************************************************/
  {
	static char temp;
	static int len,i,i1,i2,i3,j;
	static int flag3=0;

	while( (int)(temp = (char)getc(in)) != EOF)
	  {
	while((int)temp != EOF AND (int)temp != 94) temp = (char)getc(in);
	fscanf(in,"%d>     %d - %d\n",&i1,&i2,&i3);
	len=(i3-i2)/2+(i3-i2)%2;
	if(fgets(pont,len+1,in) != NULL AND i1>5 AND i3>5)
	  {
		pont[len+1]='\0';
		len = strlen(pont);
/*****************************************************************************
* Recherche des ponts : ****, **-**, *-**, **-*, **--** (ou les etoiles      *
* peuvent etre remplaces par des plus).  C.-a-d. avec l'option G-U           *
*****************************************************************************/
		for(i=0;i<len-2;i++)
		  {
		if((pont[i]=='+' OR pont[i]=='*')AND pont[i+1]!='.' AND pont[i+2] != '.' )
		  {
			if(pont[i+1]=='-' AND (pont[i+2] != '*' AND pont[i+2] != '+') AND (pont[i+3] != '*' AND pont[i+3] != '+') );
			  else
			if(pont[i+2]=='-' AND pont[i+3]=='-'AND (pont[i+4]!='*' AND pont[i+4]!='+') AND (pont[i+5]!='*' AND pont[i+5]!='+'));
			  else
				if(pont[i+3]=='-' AND pont[i+4]!='*' AND pont[i+4]!='+');
				  else
				if((pont[i]!='*' OR pont[i+1]!='*' OR pont[i+2]!='*') AND (pont[i+3]!='-' AND pont[i+3]!='+' AND pont[i+3]!='*'));
				  else
					{
					  j = i+4;
					  while(pont[j]=='+' OR pont[j]=='*')
					{
					  if(i2+j == i3-j+1)
						{
						  flag3=1;
						  break;
						}
					  j++;
					}
					  if(flag3==0)
					{
					  fprintf(out,"%d          %d   %d - %d   %d\n",i1,i,i2,i3,len);
					  fprintf(out,"%s\n",pont);
					}
					}
			break;
		  }
		  }
/*****************************************************************************
* fin de la boucle for i                                                     *
*****************************************************************************/
	  }
	  }
/*****************************************************************************
* fin de l'ensemble avec l'option G-U.                                       *
*****************************************************************************/
  }


void matrice()
/*****************************************************************************
* fonction qui gere la matrice de ponts                                      *
*****************************************************************************/
  {
	char *pont;
	char *mat;
	
	pont = (char * ) malloc(longue+2);
	if (pont==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	mat = (char *)malloc(longue+2);
	if (mat==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	if((in=fopen(fq,"r"))==NULL)
	  {
	printf("Error opening file %s for reading - end\n",fg);
	exit(1);
	  }
	if((out2 = fopen(fm,"a")) == NULL)
	  {
	printf("The file %s could not be opened.\n",fm);
	exit(3);
	  }
	mat_loop(pont,mat);
	free(pont);
	free(mat);
	fclose(in);
	fclose(out2);
  
  }


void mat_loop(char *pont, char *mat)
  {
	int i,i1,i3,i2,len;
	int k,j;

	while(fscanf(in,"%d          %d   %d - %d   %d\n",&i1,&i,&i2,&i3,&len)!=EOF)
	  {
	fscanf(in,"%s\n",pont);
	pont[len+1]='\0';
	for (j = 1; j <= longue; j++) mat[j]='0';
/*****************************************************************************
* met des zeros sur toute la longueur de la chaine                           *
*****************************************************************************/
	for (j = i; j < i+3; j++)
/*****************************************************************************
* met des 0 1 2 3  aux endroits ou il y a des bons ponts de quatre unites et    *
* plus.                                                                      *
*****************************************************************************/
	  switch(pont[j])
		{
		  case '+':
		mat[i2 + j] = '2';
		mat[i3 - j] = '2';
		break;
		  case '*':
		mat[i2 + j] = '3';
		mat[i3 - j] = '3';
		break;
		  case '-':
		mat[i2 + j] = '1';
		mat[i3 - j] = '1';
		break;
		}
	while (pont[j]=='+' OR pont[j] == '*')
	  {
		switch(pont[j])
		  {
		case '+':
		  mat[i2+j] = '2';
		  mat[i3-j] = '2';
		  break;
		case '*':
		  mat[i2 + j] = '3';
		  mat[i3 - j] = '3';
		  break;
		  }
		j++;
	  }
	if(rep=='y' OR rep=='Y')
	  {
		if(pont[j]=='-'AND (pont[j+1]=='+' OR pont[j+1]=='*'))
		  {
		mat[i2+j] = '1';
		mat[i3-j] = '1';
		j++;
		  }
		while (pont[j]=='+' OR pont[j] == '*')
		  {
		switch(pont[j])
		  {
			case '+':
			  mat[i2+j] = '2';
			  mat[i3-j] = '2';
			  break;
			case '*':
			  mat[i2 + j] = '3';
			  mat[i3 - j] = '3';
			  break;
		  }
		j++;
		  }
	  }
	  else
		{
		  if(pont[j]=='-' AND (pont[j+1]=='+' OR pont[j+1]=='*') AND (pont[j+2]=='+' OR pont[j+2]=='*'))
		{
		  mat[i2+j] = '1';
		  mat[i3-j] = '1';
		  j++;
		}
		  while (pont[j]=='+' OR pont[j] == '*')
		{
		  switch(pont[j])
			{
			  case '+':
			mat[i2+j] = '2';
			mat[i3-j] = '2';
			break;
			  case '*':
			mat[i2+j]='3';
			mat[i3-j]='3';
			break;
			}
		  j++;
		}
		}
	k = j;
	flag2=boucle_trop_petite(mat);
	if(flag2==0)
	  { 
		for (j=1;j <= longue;j++) fprintf(out2,"%c",mat[j]);
		fprintf(out2,"\n");
		dim++;
	  }
	if ((rep =='N' OR rep =='n')) mat1(pont,mat,k);
	if ((rep =='Y' OR rep =='y')) mat2(pont,mat,k);
	  }
  }
/*****************************************************************************
* Fin de la fonction mat_loop()                                              *
*****************************************************************************/


void mat1(char *pont,char *mat,int i)
/*****************************************************************************
* Fonction qui verifie s'il y a deux bons ponts ou plus sur la meme chaine   *
* de caracteres(.+*-).  Pour le choix  sans l'option G-U.                    *
*****************************************************************************/
  {
	int k;

	for (k=i; (unsigned int)k < strlen(pont)-2; k++)
	  {
	if((pont[k]=='+' OR pont[k]=='*')AND(pont[k+1]=='*'OR pont[k+1]=='+')AND pont[k+2] != '.' )
	  {
		if(pont[k+2]=='-' AND ((pont[k+3]!='*'AND pont[k+3]!='+') OR (pont[k+4]!='+' AND pont[k+4]!='*')));
		  else
		if((pont[k]!='*' OR pont[k+1]!='*' OR pont[k+2]!='*') AND (pont[k+3]!='-' AND pont[k+3]!='+' AND pont[k+3]!='*'));
		  else
			if((pont[k]!='*' OR pont[k+1]!='*' OR pont[k+2]!='*') AND pont[k+3]=='-' AND pont[k+4]!='+' AND pont[k+4]!='*'AND pont[k+5]!='+' AND pont[k+5]!='*');
			  else
			{
			  mat_loop(pont,mat);
/*****************************************************************************
* Si un deuxieme pont ou plus est trouve  on ajoute une ligne a la matrice   *
* fm avec des 0 1 2 3 a l'endroit indique                                        *
*****************************************************************************/
			  break;
			}
	  }
/*****************************************************************************
* Fin du premier if                                                          *
*****************************************************************************/
	  }
/*****************************************************************************
* Fin du for k                                                               *
*****************************************************************************/
  }
/*****************************************************************************
* Fin de la fonction mat1()                                                  *
*****************************************************************************/


void mat2(char *pont,char *mat,int k)
/*****************************************************************************
* Fonction semblable a mat1() mais pour le cas G-U.                          *
*****************************************************************************/
  {
	int i;

	for (i = k; (unsigned int)i < (strlen(pont)-2); i++)
	  {
	if((pont[i]=='+' OR pont[i]=='*')AND pont[i+1]!='.' AND pont[i+2] != '.' )
	  {
		if(pont[i+1]=='-' AND (pont[i+2] != '*' AND pont[i+2] != '+') AND (pont[i+3] != '*' AND pont[i+3] != '+') );
		  else
		if(pont[i+2]=='-' AND pont[i+3]=='-'AND (pont[i+4]!='*' AND pont[i+4]!='+') AND (pont[i+5]!='*' AND pont[i+5]!='+'));
		  else
			if(pont[i+3]=='-' AND pont[i+4]!='*' AND pont[i+4]!='+');
			  else
			if((pont[i]!='*' OR pont[i+1]!='*' OR pont[i+2]!='*') AND (pont[i+3]!='-' AND pont[i+3]!='+' AND pont[i+3]!='*'));
			  else
				{
				  mat_loop(pont,mat);
				  break;
				}
	  }
/*****************************************************************************
* Fin du premier if                                                          *
*****************************************************************************/
	  }
/*****************************************************************************
* Fin du for i                                                               *
*****************************************************************************/
  }
/*****************************************************************************
* Fin de la fonction mat2()                                                  *
*****************************************************************************/


int boucle_trop_petite(char *mat)
/*****************************************************************************
* Fonction que elimine la possibilite d'avoir des boucles infereures a la    *
* taille minimum                                                             *
*****************************************************************************/
  {
	static int c[5];
	int k1=1;
	int j1,compteur=0;

	for(j1=0;j1<2;j1++)
	  {
	while(k1<longue AND mat[k1]=='0') k1++;
	c[++compteur]=k1;
	while(k1<longue AND mat[k1]!='0') k1++;
	c[++compteur]=k1-1;
	  }
	if(c[3]-c[2] < petite-1 OR c[3]==longue OR c[2]==longue ) return(1);
	  else return(0);
  }
/*****************************************************************************
* Fin de la fonction boucle_trop_petite()                                    *
*****************************************************************************/


void mat_comp()
/*****************************************************************************
* Fonction qui cree la matrice de compatibilite                              *
*****************************************************************************/
  {
	char *temp,**mc;
	char **matcomp;
	int ii,jj,kk;

	temp = (char *)malloc(longue +2);
	if(temp==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	NEWTABLE(mc,(dim+1));
	if(mc==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	NEWTABLE(matcomp,(dim+1));
	if(matcomp==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	if ((in = fopen(fm,"r")) != NULL)
	  {
	for (ii=0;ii < dim AND (fgets(temp,longue+2,in))!=NULL;ii++)
	  {
		mc[ii] = (char *)malloc(longue+2);
		if(mc[ii]== NULL)
		  {
		printf("Allocation error - end");
		exit(2);
		  }
		strcpy(mc[ii],temp);
	  }
	fclose(in);
	  }
	  else
	{
	  printf("The file %s could not be opened.\n",fm);
	  exit(1);
	}
	printf("\n");
	for(ii=0;ii<dim;ii++)
	  {
	matcomp[ii] = (char *)malloc(dim+2);
	if (matcomp[ii]==NULL)
	  {
		printf("Allocation error - end");
		exit(2);
	  }
	for(jj=0;jj<dim;jj++)  matcomp[ii][jj] ='0';
	matcomp[ii][jj]='\0';
	  }
	for(ii=0; ii<dim; ii++)
	  for(jj = ii + 1; jj < dim; jj++)
	for(kk=0; kk < longue; kk++)
	  {
		if(mc[ii][kk]!='0' AND mc[jj][kk]!='0')
		  {
		matcomp[ii][jj]='0';
		matcomp[jj][ii]='0';
		kk = longue;
		  }
		  else
		{
		  matcomp[ii][jj] = '1';
		  matcomp[jj][ii] = '1';
		  if(rep2=='N' OR rep2=='n') croise(mc,matcomp,ii,jj);
		}
	  }
/*****************************************************************************
* Fin du for k                                                               *
*****************************************************************************/
/*****************************************************************************
* Fin du for j                                                               *
*****************************************************************************/
/*****************************************************************************
* Fin du for i                                                               *
*****************************************************************************/
	if ((out = fopen(fc,"a"))!= NULL)
	  {
	for (ii=0; ii < dim; ii++) fprintf(out,"%s \n",matcomp[ii]);
	fclose(out);
	  }
	  else
	{
	  printf("The file %s could not be opened.\n",fc);
	  exit(3);
	}
	free(temp);
	for(ii=0;ii<dim;ii++) free(mc[ii]);
	free(mc);
	kgraphe(matcomp);
	for(ii = 0; ii < dim; ii++) free(matcomp[ii]);
	free(matcomp);
  }
/*****************************************************************************
* Fin de la fonction mat_comp()                                              *
*****************************************************************************/


void croise(char **mc,char **matcomp,int i,int j)
/*****************************************************************************
* Fonction eliminant les cas croises                                         *
*****************************************************************************/
  {
	static int p[5],c[5];
	int k=0;
	int l,compteur=0;
	int compteur2=0;

	for(l=0; l<2; l++)
	  {
	while(k < longue AND mc[i][k]=='0')k++;
	p[++compteur]=k;
	while(k < longue AND mc[i][k]!='0')k++;
	c[++compteur2]=k-1;
	  }
/*****************************************************************************
* Fin du for l                                                               *
*****************************************************************************/
	k=0;
	for(l=0; l<2; l++)
	  {
	while(k < longue AND mc[j][k]=='0')k++;
	p[++compteur]=k;
	while(k < longue AND mc[j][k]!='0')k++;
	c[++compteur2]=k-1;
	  }
/*****************************************************************************
* Fin du for l                                                               *
*****************************************************************************/
	if(p[4]<p[1] OR p[2]<p[3] OR (p[1]<p[3] AND p[2]>p[4]) OR (p[3]<p[1] AND p[4]>p[2]))
	  {
	if(c[1]-p[1]==2  AND c[3]-p[3]==2) 
	  {
		matcomp[i][j] = '0';
		matcomp[j][i] = '0';
	  }  
	  }
	  else
	{
	  matcomp[i][j] = '0';
	  matcomp[j][i] = '0';
	}  
  }
/*****************************************************************************
* Fin de la fonction croisee()                                               *
*****************************************************************************/


void kgraphe(char **mcb)
/*****************************************************************************
* Fonction qui retrouve les k4 et k5 graphes a partir de la matrice de       *
* compatibilite du fichier fc...                                             *
*****************************************************************************/
  {
	int i,j,k,l,m;
	int s3,s4,s5;

	if ((out=fopen(fk,"a")) != NULL)
	  {
	for(i=0; i < dim;i++)
	  for(j=i+1; j < dim; j++)
		for(k=j+1; k < dim; k++)
		  {
		s3=0;
		s3=mcb[i][j] + mcb[i][k] + mcb[j][k] - 144;
		if(s3==3)
		  for(l=k+1; l < dim; l++)
			{
			  s4=0;
			  s4=s3 + mcb[i][l] + mcb[j][l] + mcb[k][l] - 144;
			  if(s4==6)
			{
			  printf("k4 :%d - %d - %d - %d - %d\n",i+1,j+1,k+1,l+1,0);
			  fprintf(out,"k4 :%d , %d , %d , %d , %d\n",i+1,j+1,k+1,l+1,0);
			  knombre++;
			  k4count++;
			  for(m=l+1; m < dim; m++)
				{
				  s5=0;
				  s5=s4 + mcb[i][m] + mcb[j][m] + mcb[k][m] + mcb[l][m] -192;
				  if(s5==10)
				{
				  printf("k5 :%d - %d - %d - %d - %d\n",i+1,j+1,k+1,l+1,m+1);
				  fprintf(out,"k5 :%d , %d , %d , %d , %d\n",i+1,j+1,k+1,l+1,m+1);
				  knombre++;
				  k5count++;
				}
				}
/*****************************************************************************
* Fin du for m                                                               *
*****************************************************************************/
			}
/*****************************************************************************
* Fin du if s4                                                               *
*****************************************************************************/
			}
/*****************************************************************************
* Fin du for l                                                               *
*****************************************************************************/
/*****************************************************************************
* Fin du if s3                                                               *
*****************************************************************************/
		  }
/*****************************************************************************
* Fin du for k                                                               *
*****************************************************************************/
/*****************************************************************************
* Fin du for j                                                               *
*****************************************************************************/
/*****************************************************************************
* Fin du for i                                                               *
*****************************************************************************/
	fclose(out);
	  }
/*****************************************************************************
* Fin du if fopen(fk)                                                        *
*****************************************************************************/
	  else
	{
	  printf("Error attempting to open file %s - end.\n",fk);
	  exit(3);
	}
  }
/*****************************************************************************
* Fin de la procedure kgraphe                                                *
*****************************************************************************/


void fl_fichier()
/*****************************************************************************
* Fonction qui cree le fichier fl...  C.a.d. le fichier contenant les        *
* numeros identifiants la molecule et les lignes corresondantes a ces        *
* numeros a partir de la matrice du fichier fm...                            *
*****************************************************************************/
  {
	int i,j,j1=0,j2=0;
	int k,i1,i2,i3,i4,i5;
	int flag=0;
	char *pont,temp;

	i1=0;i2=0;i3=0;i4=0;i5=0;
	pont = (char *)malloc(longue +2);
	if (pont==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	if ((in = fopen(fk,"r")) == NULL)
	  {
	printf("Error opening file %s for reading - end\n",fk);
	exit(1);
	  }
	if((out = fopen(fl,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",fl);
	exit(3);
	  }
	if((in1 = fopen(fm,"r")) == NULL)
	  {
	printf("Error opening file %s for reading - end\n",fm);
	exit(1);
	  }
	for(i=1; i <= knombre; i++)
	  {
	temp=(char)fgetc(in);
	while((int)temp!= EOF AND (int)temp!=58) temp=(char)fgetc(in);
	fscanf(in,"%d , %d , %d , %d , %d\n",&i1,&i2,&i3,&i4,&i5);
	fprintf(out,"%d , %d , %d , %d , %d\n",i1,i2,i3,i4,i5);
	for(k=1; k <= 5; k++)
	  {
		switch(k)
		  {
		case 1:
		  j1=0;
		  j2=i1;
		  break;
		case 2:
		  j1=i1;
		  j2=i2;
		  break;
		case 3:
		  j1=i2;
		  j2=i3;
		  break;
		case 4:
		  j1=i3;
		  j2=i4;
		  break;
		case 5:
		  j1=i4;
		  j2=i5;
		  break;
		  }
		for(j=j1; j < j2; j++)
		  {
		flag=1;
		if(fgets(pont,longue+2,in1) == NULL)
		  {
			printf("Error reading file %s - end\n",fm);
			exit(4);
		  }
		  }
		if(flag==1) fprintf(out,"%s",pont);
		flag=0;
	  }
/*****************************************************************************
* Fin du for k                                                               *
*****************************************************************************/
	rewind(in1);
	  }
/*****************************************************************************
* Fin du for i                                                               *
*****************************************************************************/
	free(pont);
	fclose(in1);
	fclose(out);
	fclose(in);
  }
/*****************************************************************************
* Fin de la fonction fl_fichier()                                            *
*****************************************************************************/


void fd_fichier()
/*****************************************************************************
* Fonction qui lit le fichier fl... et cree le fichier fd...                 *
*****************************************************************************/
  {
	int i1,i2,i3,i4,i5,j,l;
	char *tmp, **ligne;

	tmp=(char *)malloc(longue+2);
	if(tmp==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	if((out=fopen(fd,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",fd);
	exit(3);
	  }
	if((in = fopen(fl,"r"))==NULL)
	  {
	printf("Error opening file %s to create file %s - fin\n",fl,fd);
	exit(1);
	  }
	for(l=1;l<=knombre;l++)
	  {
	fscanf(in,"%d , %d , %d , %d , %d\n",&i1,&i2,&i3,&i4,&i5);
	if(i5==0) nboucle=4; 
	  else nboucle=5;
	NEWTABLE(ligne,(nboucle+1));
	if (ligne==NULL)
	  {
		printf("Allocation error - end");
		exit(2);
	  }
	for(j=1; j <= nboucle AND (fgets(tmp,longue+2,in))!=NULL;j++)
	  {
		ligne[j] = (char *)malloc(longue+2);
		if(ligne[j] == NULL)
		  {
		printf("Allocation error - end");
		exit(2);
		  }
		strcpy(ligne[j],tmp);
	  }
	fprintf(out,"> %d , %d , %d , %d , %d\n",i1,i2,i3,i4,i5);
	marque_pont(ligne);
	plus_petit();
	couplage();
	calcul_pont_boucle(ligne);
	for(j=1;j<=nboucle;free(ligne[j++]));
	free(ligne);
	if (rep2=='N' OR rep2=='n') dessin(i1,i2,i3,i4,i5);
	ncouple=0;
	  }
	free(tmp);
	fclose(in);
	fclose(out);
  }
/*****************************************************************************
* Fin de la fonction fd_fichier().                                           *
*****************************************************************************/


void marque_pont(char **lp)
/*****************************************************************************
* Fonction qui marque le debut et la fin de chaque pont                      *
*****************************************************************************/
  {
	int i,j,k,temp,compteur=0;

	k=0;
	for(i=1;i<=nboucle;i++)
	  {
	for(j=1;j<=2;j++)
	  {
		while(k<longue AND lp[i][k]=='0')k++;
		if(k!=longue)
		  {
		pl[++compteur]=k+1;
		while(k<longue AND lp[i][k]!='0')k++;
		pl[++compteur]=k;
		  }
	  }
	if((compteur%4)!=0)
	  {
		temp=pl[compteur];
		pl[compteur]=(temp-pl[compteur-1]+1)/2 +pl[compteur-1]-1;
		pl[++compteur]=pl[compteur]+1;
		pl[++compteur]=temp;
	  }
	k=0;
	  }
  }
/*****************************************************************************
* Fin de la fonction marque_pont().                                          *
*****************************************************************************/


void plus_petit()
/*****************************************************************************
* Fonction qui met en ordre croissante dans la variable pp[], les ponts      *
* premiers sur chaque ligne impliquee dans la molecule, afin de faire le     *
* couplage, le calcul des ponts et boucles, etc...                           *
*****************************************************************************/
  {
	int c[21];
	int indice=0;
	int temp,i,j;

	for(i=1; i<=4*nboucle;i++) c[i]=pl[i];
	temp=c[1];
	for(j=1; j<=nboucle;j++)
	  {
	for(i=1;i<=nboucle*4; i=i+4) if(temp>c[i]) temp=c[i];
	pp[++indice]=temp;
	for(i=1; i<=4*nboucle;i=i+4)
	  {
		if(pp[indice]==c[i])
		  {
		c[i]=longue+1;
		temp=c[i];
		break;
		  }
	  }
	  }
	pp[++indice]=longue;
  }
/*****************************************************************************
* Fin de la fonction plus_petit()                                            *
*****************************************************************************/


void couplage()
/*****************************************************************************
* Fonction qui relie les positions des liens deux a deux                     *
*****************************************************************************/
  {
	int i,j,k,l,m,n,temp,temp1,temp2;
	int indice=1;
	int flag=0;
	int compteur=0;

	i=0; j=0; k=0; l=0; m=0; n=0;
	while(++j<pp[indice])
	  {
	ecrire_fd(j,0);
	compteur++;
	  }
	while(compteur<longue)
	  {
	while(compteur>pp[indice]) indice++;
	for(l=1;l<=4*nboucle;l=l+2)
	  {
		if(pl[l]==pp[indice] AND pl[l]==compteur+1)
		  {
		flag=1;
		break;
		  }
		  else
		if(pl[l]==compteur+1)
		  {
			flag=2;
			break;
		  }
	  }
	if(l==4*nboucle+1) l=m;
	for(j=pl[l],k=pl[l+3]; j<=pl[l+1] AND flag==1;j++,k--)
	  {
		ecrire_fd(j,k);
		compteur++;
	  }
	for(j=pl[l],k=pl[l-1];j<=pl[l+1] AND flag==2;j++,k--)
	  {
		ecrire_fd(j,k);
		compteur++;
	  }
	for(m=1;m<=4*nboucle;m=m+2)
	  if(pl[m]>pl[l])
		{
		  for(n=1;n<=4*nboucle;n=n+2)
		if(pl[n]<pl[m] AND pl[n]>pl[l]) m=n;
		  break;
		}
	if(m==4*nboucle+1)
	  {
		m=l+2;
		temp=longue+1;
	  }
	  else temp=pl[m];
	if(l>=4*nboucle-1) temp1=longue+1;
	  else
		if (((l+1)%4)==0) temp1=temp;
		  else temp1=pl[l+2];
	if(indice==nboucle+1) temp2=pp[indice]+1;
	  else temp2=pp[indice+1];
	for(k=compteur+1;k<temp2 AND k<temp1 AND k<temp;k++)
	  {
		ecrire_fd(k,0);
		compteur++;
	  }
	if(k==pl[l+2]) for(j=pl[l+2],k=pl[l+1];j<=pl[l+3];j++,k--)
	  {
		ecrire_fd(j,k);
		compteur++;
	  }
	flag=0;
	  }
  }
/*****************************************************************************
* Fin de la fonction couplage()                                              *
*****************************************************************************/


void ecrire_fd(int l,int m)
/*****************************************************************************
* Fonction qui ecrit le couplage sur fd...                                   *
*****************************************************************************/
  {
	if(ncouple==10)
	  {
	fprintf(out,"\n");
	ncouple=0;
	  }
	fprintf(out,"%d.%d , ",l,m);
	ncouple++;
  }
/*****************************************************************************
* Fin de la fonction ecrire_fd().                                            *
*****************************************************************************/


void calcul_pont_boucle(char **lp)
/*****************************************************************************
* Fonction qui calcule le nombre de ponts dans la molecule et le nombre de   *
* bases impliquees dans chaque boucle.                                       *
*****************************************************************************/
  {
	int pont=0;
	int maille=0;
	char temp,b_type=' ';
	float enerpont=0.0F;
	int i,j,k,l,m,x;

	for(i=1,j=1;j<=nboucle AND i<=4*nboucle;i+=4,j++)
	  {
	pont=pont+(pl[i+1]-pl[i]+1);
	k=0;
	while(k<=longue AND lp[j][k]=='0')k++;
	while(k<=longue AND lp[j][k]!='0')
	  {
		temp=lp[j][k];
		switch(temp)
		  {
		case '1':
		  enerpont += 6.0F;
		  break;
		case '2':
		  enerpont += 11.7F;
		  break;
		case '3':
		  enerpont += 18.0F;
		  break;
		  }
		k++;
	  }
	  }
	fprintf(out,"\n  BRIDGES = %d , Bridge energy = %0.1f \n",pont,enerpont);
	if (rep2=='n' OR rep2=='N')
	  {
	fprintf(out,"  LOOPS = ");
	for(k=nboucle;k>=1;k--)
	  {
		for(l=1;l<=4*nboucle;l=l+2) if(pp[k]==pl[l]) break;
		if(k==nboucle)
		  {
		b_type='H';
		maille=pl[l+2]-pl[l+1]+1;
		  }
		  else
		{
		  x=retrouve_pont_interieur(l);
		  m=l;
		  if(x==l+2)
			{
			  b_type='H';
			  maille=pl[l+2]-pl[l+1]+1;
			}
			else
			  {
			maille=pl[x]-pl[m+1]+1;
			do
			  {
				for(j=1;j<=nboucle;j++) if(pp[j]==pl[x])break;
				if (j==nboucle+1)
				  {
				m=x;
				x=retrouve_pont_interieur(x);
				maille=maille +pl[x]-pl[m+1]+1;
				if((pl[m-2]-pl[x-1]==1) OR(pl[x]-pl[m+1]==1)) b_type='B'; else b_type='I';
				  }
				  else x=x+2;
			  }
			while(x!=l+2);
			  }
		}
		fprintf(out,"%d %c  ,  ",maille,b_type);
		maille=0;
	  }
	  }
/*****************************************************************************
* Fin de la condition non-croisee                                            *
*****************************************************************************/
	fprintf(out,"\n\n");
  }
/*****************************************************************************
* Fin de la fonction calcul_pont_boucle()                                    *
*****************************************************************************/


int retrouve_pont_interieur(int i)
/*****************************************************************************
* Fonction qui verifie s'il y a un pont a l'interieur du pont actuelement    *
* etudie par la fonction calcul_pont_boucle()                                *
*****************************************************************************/
  {
	int m,n;

	for(m=1;m<=4*nboucle;m=m+2)
	  if (pl[m]>pl[i])
	{
	  for(n=1;n<=4*nboucle;n=n+2)
		if(pl[n]<pl[m] AND pl[n]>pl[i]) m=n;
	  break;
	}
	return(m);
  }
/*****************************************************************************
* Fin de la fonction retrouve_pont_interieur()                               *
*****************************************************************************/


void dessin(int j1,int j2,int j3,int j4,int j5)
/*****************************************************************************
* Fonction qui esquisee la forme de certaines molecules telle que            *
* les trefles, lineaires, etc..                                              *
*****************************************************************************/
  {
	int i,j,l,pg=0;
	int Tailtype=0,indice=0;
	int cond[7];

	if((out1=fopen(fe,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",fe);
	exit(3);
	  }
	fprintf(out1,"\n %d , %d , %d , %d , %d\n\n",j1,j2,j3,j4,j5);
	for(j=1;j<=4*nboucle;j++) if(pg<pl[j]) pg=pl[j];
	if(nboucle==4)
	  {
	if(pp[1] > 1 AND pg >= longue)
	  {
		fprintf(out1,"               |\n               |\n               ||\n               ||\n");
		Tailtype=1;
	  }
	if(pg < longue AND pp[1] <= 1)
	  {
		fprintf(out1,"                |\n                |\n               ||\n               ||\n");
		Tailtype=2;
	  }
	if(pg < longue AND pp[1] > 1)
	  {
		fprintf(out1,"             \\    /\n              \\  /\n               ||\n               ||\n");
		Tailtype=3;
	  }
	if(pp[1] <= 1 AND pg >= longue)
	  {
		fprintf(out1,"               ||\n               ||\n");
		Tailtype=4;
	  }
	for(i=1;i<=nboucle;i++)
	  {
		for(l=1;l<=4*nboucle;l=l+4) 
		  if(pl[l]==pp[i])
		{
		  cond[++indice]=pl[l+2];
		  break;
		}   
	  }
	if(cond[4]<cond[3] AND cond[3]<cond[2] AND cond[2]<cond[1])
	  {
		fprintf(out1," O====O====O====O\n\n");
		counts[Tailtype][1]++;
	  }
	  else
		if(cond[2]<cond[3] AND cond[3]<cond[4] AND cond[4]<cond[1] AND cond[2]<pp[3] AND cond[3]<pp[4])
		  {
		fprintf(out1,"          O====O====O\n");
		fprintf(out1,"               ||\n               ||\n               O\n\n");
		counts[Tailtype][2]++;
		  }
		  else
		if(cond[3]<cond[4] AND cond[4]<cond[2] AND cond[2]<cond[1] AND cond[3]<pp[4])
		  {
			fprintf(out1,"     O====O====O\n");
			fprintf(out1,"         ||\n         ||\n         O\n\n");
			counts[Tailtype][3]++;
		  }
		  else
			if(cond[3]<cond[2] AND cond[2]<cond[4] AND cond[4]<cond[1] AND cond[2]<pp[4])
			  {
			fprintf(out1,"     O====O====O====O\n\n");
			counts[Tailtype][4]++;
			  }
			  else
			if(cond[2]<cond[4] AND cond[4]<cond[3] AND cond[3]<cond[1] AND cond[2]<pp[3])
			  {
				fprintf(out1,"          O====O====O====O\n\n");
				counts[Tailtype][5]++;
			  }
			  else 
				{
				  fprintf(out1,"Unrecognized condition \n\n");
				  counts[Tailtype][18]++;
				}
	  }
	  else if(nboucle==5)
	{
	  if(pp[1] > 1 AND pg >= longue)
		{
		  fprintf(out1,"                    |\n                    |\n                    ||\n                    ||\n");
		  Tailtype=1;
		}
	  if(pg < longue AND pp[1] <= 1) 
		{
		  fprintf(out1,"                     |\n                     |\n                    ||\n                    ||\n");
		  Tailtype=2;
		}
	  if(pg < longue AND pp[1] > 1) 
		{
		  fprintf(out1,"                  \\    /\n                   \\  /\n                    ||\n                    ||\n");
		  Tailtype=3;
		}
	  if(pp[1] <= 1 AND pg >= longue) 
		{
		  fprintf(out1,"                    ||\n                    ||\n");
		  Tailtype=4;
		}
	  for(i=1;i<=nboucle;i++)
		{
		  for(l=1;l<=4*nboucle;l=l+4) 
		if(pl[l]==pp[i])
		  {
			cond[++indice]=pl[l+2];
			break;
		  }   
		}
	  if(cond[5]<cond[4] AND cond[4]<cond[3] AND cond[3]<cond[2] AND cond[2]<cond[1])
		{
		  fprintf(out1," O====O====O====O====O\n\n");
		  counts[Tailtype][6]++;
		}
		else
		  if(cond[4]<cond[3] AND cond[3]<cond[2] AND cond[2]<cond[5] AND cond[5]<cond[1] AND cond[2]<pp[5])
		{
		  fprintf(out1,"      O====O====O====O====O\n\n");
		  counts[Tailtype][7]++;
		}
		else
		  if(cond[2]<cond[5] AND cond[5]<cond[4] AND cond[4]<cond[3] AND cond[3]<cond[1] AND cond[2]<pp[3])
			{
			  fprintf(out1,"                O====O====O====O====O\n\n");
			  counts[Tailtype][8]++;
			}
			else
			  if(cond[3]<cond[2] AND cond[2]<cond[5] AND cond[5]<cond[4] AND cond[4]<cond[1] AND cond[2]<pp[4])
			{
			  fprintf(out1,"           O====O====O====O====O\n\n");
			  counts[Tailtype][9]++;
			}
			else
			  if(cond[4]<cond[3] AND cond[3]<cond[5] AND cond[5]<cond[2] AND cond[2]<cond[1] AND cond[3]<pp[5])
				{
				  fprintf(out1,"      O====O====O====O\n");
				  fprintf(out1,"                ||\n                ||\n                O\n\n");
				  counts[Tailtype][10]++;
				}
				else
				  if(cond[3]<cond[5] AND cond[5]<cond[4] AND cond[4]<cond[2] AND cond[2]<cond[1] AND cond[3]<pp[4])
				{
				  fprintf(out1,"                     O====O====O====O\n");
				  fprintf(out1,"                          ||\n                          ||\n                          O\n\n");
				  counts[Tailtype][11]++;
				}
				else
				  if(cond[3]<cond[2] AND cond[2]<cond[4] AND cond[4]<cond[5] AND cond[5]<cond[1] AND cond[2]<pp[4] AND cond[4]<pp[5])
					{
					  fprintf(out1,"           O====O====O====O\n");
					  fprintf(out1,"                     ||\n                     ||\n                     O\n\n");
					  counts[Tailtype][12]++;
					}    
					else
					  if(cond[2]<cond[4] AND cond[4]<cond[5] AND cond[5]<cond[3] AND cond[3]<cond[1] AND cond[2]<pp[3] AND cond[4]<pp[5])
					{
					  fprintf(out1,"                O====O====O====O\n");
					  fprintf(out1,"                          ||\n                          ||\n                          O\n\n");
					  counts[Tailtype][13]++;
					}
					else
					  if(cond[3]<cond[4] AND cond[4]<cond[2] AND cond[2]<cond[5] AND cond[5]<cond[1] AND cond[3]<pp[4] AND cond[2]<pp[5])
						{
						  fprintf(out1,"           O====O====O====O\n");
						  fprintf(out1,"                ||\n                ||\n                O\n\n");
						  counts[Tailtype][14]++;
						}
						else
						  if(cond[4]<cond[5] AND cond[5]<cond[3] AND cond[3]<cond[2] AND cond[2]<cond[1] AND cond[4]<pp[5])
						{
						  fprintf(out1,"                     O====O====O====O\n");
						  fprintf(out1,"                               ||\n                               ||\n                               O\n\n");
						  counts[Tailtype][15]++;
						}
						else
						  if(cond[3]<cond[4] AND cond[4]<cond[5] AND cond[5]<cond[2] AND cond[2]<cond[1] AND cond[3]<pp[4] AND cond[4]<pp[5])
							{
							  fprintf(out1,"           O====O====O\n");
							  fprintf(out1,"               //\\\\\n              //  \\\\\n             O      O\n\n");
							  counts[Tailtype][16]++;
							}
							else
							  if(cond[2]<cond[3] AND cond[3]<cond[4] AND cond[4]<cond[5] AND cond[5]<cond[1] AND cond[2]<pp[3] AND cond[3]<pp[4] AND cond[4]<pp[5])
							{
							  fprintf(out1,"                O====O====O\n");
							  fprintf(out1,"                    //\\\\\n                   //  \\\\\n                  O      O\n\n");
							  counts[Tailtype][17]++;
							}
							else 
							  {
								fprintf(out1,"Unrecognized condition\n\n");
								counts[Tailtype][18]++;
							  }
	}
	else 
	  if(nboucle<4) fprintf(out1,"Cases with less than 4 loops were not treated.\n\n");
		else
		  if(nboucle>5) fprintf(out1,"Cases with more than 5 loops were not treated.\n\n");
		else fprintf(out1,"Unrecognized condition\n\n");
	fclose(out1);
  }
/*****************************************************************************
* Fin de la fonction dessin()                                                *
*****************************************************************************/


void ft_fichier()
/*****************************************************************************
* Fonction qui creer le fichier ft... C.a.d le tableau contenant les numeros *
* identifiant chaque molecule, leur nombre de ponts, leurs nombres de bases  *
* par boucle et l'energie impliquee.                                         *
*****************************************************************************/
  {
	int i1,i2,i3,i4,i5,j1,j2,j3,j4,j5;
	char h1,h2,h3,h4,h5;
	int pont;
	float energie,e_pont;
	char temp;
	static char ligne[]="---------------------------------------------------------------------------------";
	static char cote[]="|                                                                               |";

	if((out=fopen(ft,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",ft);
	exit(3);
	  }
	if((in=fopen(fd,"r"))==NULL)
	  {
	printf("Error opening file %s - end\n",fd);
	exit(1);
	  }
	fprintf(out,"%s\n",ligne);
	fprintf(out,"%s\n",cote);
	fprintf(out,"|                                         RNA%s                                |\n",name);
	if(rep=='n' OR rep=='N')
	  fprintf(out,"|                                   WITHOUT NON-WATSON-CRICK                    |\n");
	  else
	fprintf(out,"|                                   WITH NON-WATSON-CRICK                       |\n");
	fprintf(out,"%s\n",cote);
	fprintf(out,"%s\n",ligne);
	fprintf(out,"|  IDENTIFICATION   |BRIDGES|          LOOPS        |          ENERGY           |\n");
	fprintf(out,"|                   |       |                       |   Loops  | Bridges| L - B |\n");
	fprintf(out,"%s\n",ligne);
	while( (int)(temp=(char)getc(in)) != EOF)
	  {
	while((int)temp!=EOF AND (int)temp!=62) temp=(char)getc(in);
	fscanf(in," %d , %d , %d , %d , %d\n",&i1,&i2,&i3,&i4,&i5);
	while((int)temp!=EOF AND (int)temp!=66) temp=(char)getc(in);
	fscanf(in,"RIDGES = %d , Bridge energy = %f\n",&pont,&e_pont);
	if(i5==0) fscanf(in,"  LOOPS = %d %c  ,  %d %c  ,  %d %c  ,  %d %c  ,\n",&j1,&h1,&j2,&h2,&j3,&h3,&j4,&h4);
	  else fscanf(in,"  LOOPS = %d %c  ,  %d %c  ,  %d %c  ,  %d %c ,  %d %c  ,\n",&j1,&h1,&j2,&h2,&j3,&h3,&j4,&h4,&j5,&h5);
	if (i5==0) 
	  {
		j5=0;
		h5=' '; 
	  }
	energie = calcul_energie(j1,h1,j2,h2,j3,h3,j4,h4,j5,h5);
	fprintf(out,"|%2d -%2d -%2d -%2d -%2d |  %d   |",i1,i2,i3,i4,i5,pont);
	if (i5==0)
	  fprintf(out,"%2d%c  %2d%c  %2d%c  %2d%c     |  %4.2f  | %3.2f |%3.2f |\n",j1,h1,j2,h2,j3,h3,j4,h4,energie,e_pont,energie-e_pont);
	  else
		fprintf(out,"%2d%c  %2d%c  %2d%c  %2d%c  %2d%c|  %4.2f  | %3.2f |%3.2f |\n",j1,h1,j2,h2,j3,h3,j4,h4,j5,h5,energie,e_pont,energie-e_pont);
	  }
	fprintf(out,"%s\n",ligne);
	fclose(in);
	fclose(out);
  }
/*****************************************************************************
* Fin de la fonction ft_fichier()                                            *
*****************************************************************************/


float calcul_energie(int j,char j1,int k,char k1,int l,char l1,int m,char m1,int n,char n1)
/*****************************************************************************
* Fonction qui calcule l'energie de la molecule transmise par la fonction    *
* ft_fichier()                                                               *
*****************************************************************************/
  {
	float x,y,z,v,w,energie;

	x=nombre_boucle(j,j1);
	y=nombre_boucle(k,k1);
	z=nombre_boucle(l,l1);
	v=nombre_boucle(m,m1);
	w=nombre_boucle(n,n1);
	energie= v + w + x + y + z;
	return(energie);
  }
/*****************************************************************************
* Fin de la fonction calcul_energie()                                        *
*****************************************************************************/


float nombre_boucle(int i,char i1)
/*****************************************************************************
* Fonction qui attribue le montant d'energie correspondante au nombre de     *
* bases impliquees dans la boucle transmise par la fonction                  *
* calcul_energie().                                                          *
*****************************************************************************/
  {
	switch(i)
	  {
	case 0:
	  return(0.0F);
	  break;
	case 2:
	  if(i1=='B') return(100.0F);
		else
		  if(i1=='I') return(10.0F);
		else return(350.0F);
	  break;
	case 3:
	  if(i1=='B') return(176.95F);
		else
		  if(i1=='I') return(34.65F);
		else return(326.8F);
	  break;
	case 4:
	  if(i1=='B') return(196.20F);
		else
		  if(i1=='I') return(65.40F);
		else return(196.08F);
	  break;
	case 5:
	  if(i1=='B') return(207.7F);
		else
		  if(i1=='I') return(84.65F);
		else return(161.52F);
	  break;
	case 6:
	  if(i1=='B') return(215.4F);
		else
		  if(i1=='I') return(100.0F);
		else return(169.2F);
	  break;
	case 7:
	  if(i1=='B') return(219.25F);
		else
		  if(i1=='I') return(111.55F);
		else return(176.88F);
	  break;
	case 8:
	  if(i1=='B') return(221.2F);
		else
		  if(i1=='I') return(115.4F);
		else return(184.56F);
	  break;
	case 9:
	  if(i1=='B') return(223.1F);
		else
		  if(i1=='I') return(117.4F);
		else return(188.4F);
	  break;
	case 10:
	  if(i1=='B') return(226.95F);
		else
		  if(i1=='I') return(119.25F);
		else return(196.08F);
	  break;
	case 11:
	  if(i1=='B') return(230.8F);
		else
		  if(i1=='I') return(123.1F);
		else return(200.0F);
	  break;
	case 12:
	  if(i1=='B') return(234.5F);
		else
		  if(i1=='I') return(125.0F);
		else return(203.8F);
	  break;
	case 13:
	  if(i1=='B') return(236.5F);
		else
		  if(i1=='I') return(126.95F);
		  else return(207.6F);
	  break;
	case 14:
	  if(i1=='B') return(238.5F);
		else
		  if(i1=='I') return(130.8F);
		else return(208.0F);
	  break;
	case 15:
	  if(i1=='B') return(238.5F);
		else
		  if(i1=='I') return(130.8F);
		else return(211.8F);
	  break;
	default:
	  if(i1=='B') return(238.5F);
		else
		  if(i1=='I') return(130.8F);
		else return(211.8F);
	  break;
	  }
  }
/*****************************************************************************
* Fin de la procedure nombre_boucle()                                        *
*****************************************************************************/


void fx_fichier()
/*****************************************************************************
* Fonction qui cree le fichier fx... C.a.d le fichier contenant le numero    *
* de lignes correspondant a l'identification des molecules suivit de la      *
* description en nombres de leurs ponts                                      *
*****************************************************************************/
  {
	int i,j,k,marque[5];
	char *temp;
	int compteur=0;

	if((in=fopen(fm,"r"))==NULL)
	  {
	printf("Error opening file %s - end\n",fm);
	exit(1);
	  }
	if((out=fopen(fx,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",fx);
	exit(3);
	  }
	fprintf(out,"                         RNA%s \n",name);
	if (rep=='N' OR rep=='n') fprintf(out,"                   WITHOUT NON-WATSON-CRICK\n");
	  else fprintf(out,"                   WITH NON-WATSON-CRICK\n\n");
	temp = (char *)malloc(longue+2);
	if(temp==NULL)
	  {
	printf("Allocation error - end");
	exit(2);
	  }
	for(i=1;i<=dim AND (fgets(temp,longue+2,in))!=NULL;i++)
	  {
	k=0;
	compteur=0;
	for(j=1;j<=2;j++)
	  {
		while(k<longue AND temp[k]=='0')k++;
		marque[++compteur]=k+1;
		while(k<longue AND temp[k]!='0')k++;
		marque[++compteur]=k;
	  }
	fprintf(out,"%2d.    ",i);
	for(j=marque[1],k=marque[4]; j<=marque[2];j++,k--)
	fprintf(out,"  %2d.%2d  ,",j,k);
	fprintf(out,"\n");
	  }
	free(temp);
	fclose(in);
	fclose(out);
  }
/*****************************************************************************
* Fin de la fonction fx_fichier()                                            *
*****************************************************************************/


void ff_fichier()
/*****************************************************************************
* Fonction qui cree le fichier ff... C.a.d le fichier contenant les          *
* diagrammes de connectivite                                                 *
*****************************************************************************/
  {
	int i9,j9,k9,i1,i2,i3,i4,i5;
	int compteur=0;
	char temp;

	if((in=fopen(fd,"r"))==NULL)
	  {
	printf("Error opening file %s - end\n",fd);
	exit(1);
	  }
	if((out=fopen(ff,"a"))==NULL)
	  {
	printf("Error opening file %s - end.",ff);
	exit(3);
	  }
	while( (int)(temp=(char)getc(in)) != EOF)
	  {
	while((int)temp!=EOF AND (int)temp!=62) temp=(char)getc(in);
	if((int)temp==EOF) break;
	fscanf(in," %d , %d , %d , %d , %d\n",&i1,&i2,&i3,&i4,&i5);
	fprintf(out,"%d , %d , %d , %d , %d\n",i1,i2,i3,i4,i5);
	for(i9=1;i9<=longue;i9++)
	  {
		fscanf(in,"%d.%d , ",&j9,&k9);
		compteur++;
		if(compteur==10)
		  {
		fscanf(in,"\n");
		compteur=0;
		  }
		if(k9==0) fprintf(out,"%2d %c\n",j9,chaine[j9-1]);
		  else
		{
		  fprintf(out,"%2d ",j9);
		  fprintf(out,"%c - ",chaine[j9-1]);
		  fprintf(out,"%c ",chaine[k9-1]);
		  fprintf(out,"%2d\n",k9);
		}
	  }
	fprintf(out,"\n");
	  }
	fclose(in);
	fclose(out);
	free(chaine);
  }
/*****************************************************************************
* Fin de la fonction ff_fichier()                                            *
*****************************************************************************/


void fp_fichier()
/*****************************************************************************
* Fonction qui cree le fichier fp... C.a.d le fichier contenant les          *
* comptes des K graphs et des structures                                     *
*****************************************************************************/
  {
	int i,j;

	if((out=fopen(fp,"a"))==NULL)
	  {
	printf("Error opening file %s - end\n",fp);
	exit(3);
	  }
	fprintf(out,"There are %d K graphs of order 4 and %d K graphs of order 5.\n\n",k4count,k5count);
	if(rep2=='N' OR rep2=='n')
	  {
	fprintf(out,"  Drawing code     Count\n------------------------\n");
	for(i=1;i<5;i++)
	  for(j=1;j<19;j++)
		fprintf(out,"       %c%c             %d\n",(char)(i+64),(char)(j+64),counts[i][j]);
	fprintf(out,"------------------------\n");
	  }
	fclose(out);
  }
/*****************************************************************************
* Fin de la fonction fp_fichier()                                            *
*****************************************************************************/


void main()
  {
	trna1();
	init_fichiers();
	invers();
	glisse();
	cherche_pont();
	matrice();
	mat_comp();
	if(knombre!=0)
	  {
	fl_fichier();
	fd_fichier();
	if(rep2=='n' OR rep2=='N') ft_fichier();
	fx_fichier();
	ff_fichier();
	fp_fichier();
	  }
	  else
	printf("There are no K graphs of order 4 or 5.\n");
	fcloseall();
	exit(0);
  }
/*****************************************************************************
*                                                                            *
*                       ********  *******  *       *                         *
*                       *            *     * *     *                         *
*                       *****        *     *   *   *                         *
*                       *            *     *     * *                         *
*                       *         *******  *       *                         *
*                                                                            *
*****************************************************************************/
