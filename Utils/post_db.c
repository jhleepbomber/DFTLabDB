#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h> 
#include <stdbool.h>

#define BUFSIZE 1000

#define N_nod_mx 140000000
#define N_atoms_mx 15000
#define N_bands_type_mx 15000
#define N_bands_mx 1000
#define N_DOS_mx 2000

#define Max_Atom_name 1000

#define MAX_STR_LEN 4000

#define Nx_ 500
#define Ny_ 500
#define Nz_ 500

#define  MAX_RhoV 3000000


char*  trim(char *s); // ���ڿ� �¿� ���� ��� ���� �Լ�
char* ltrim(char *s); // ���ڿ� ���� ���� ���� �Լ�
char* rtrim(char* s); // ���ڿ� ���� ���� ���� �Լ�



double func[N_nod_mx];

int  N_P_Atom[N_atoms_mx];
double  P_Atom[N_atoms_mx][4];
	char   ID_Atom[N_atoms_mx][Max_Atom_name]={};	
	char   Link[N_atoms_mx][Max_Atom_name]={};
	char   ID_cAtom[N_atoms_mx][Max_Atom_name]={};


	
	int OverLap_Atom[N_atoms_mx];

	int N_orbit_per_atom[N_atoms_mx];

double Band_kE[N_bands_type_mx][N_bands_mx][2];
int    Band_spin[N_bands_type_mx][N_bands_mx];

double DOS[N_DOS_mx][3];
double E_pDOS[N_DOS_mx];

double pDOS_T[N_DOS_mx];
double pDOS_E_T[N_DOS_mx];

char args1[230]="";
char SystemLabel[230]="";

int N_nod, N_meshx;	

int Num_atom(char* );	

int Nx, Ny, Nz, Nzy;

int N_atoms=0, N_V_atoms=0, N_Links=0;


double dx, dy, dz;

int xi=0, yi=0, zi=0;

int Ncei[256];
int ei[256][4][3];
int iei[13][2][3];

int orderF(int, int, int);

//int ei[1][5]={{0, 0, 0, 0, 0}};

	double xmx=-1e+100, ymx=-1e+100, zmx=-1e+100;
    double xmn= 1e+100, ymn= 1e+100, zmn= 1e+100;

	double func_mx=-1e+100, func2_mx=-1e+100;
	double func_mn= 1e+100, func2_mn= 1e+100;



	int Is_same;	

	int c_i=0;
	int cL_i=0;

	int md;

	int N_t=1; 	
	int ti;

	int F_DOS=0;


	char On_save_R[10]="";
	char On_save_P[10]="";

double absol(double value);

void Save_Bonds(char*);

void VR_info(char*);
int Save_V(char*), Save_Rho(char*);
void Save_S_V_R(char*, char*);


void Eliminate(char *str, char ch);

char *replaceAll(char *s, const char *olds, const char *news);

int point_count(char *);

int main(int argc, char **args)
{
	
	FILE *fr, *fw, *fr_data, *fw_data;
	FILE *fwT;
//	int state0;

	//FILE *savedata0, *savedata1,*savedata2,*savedata3,*savedata4,*savedata5,*savedata6;
	

	FILE *read_data0, *read_data2;
	

	char filename[230]="";
	char filename2[230]="";
	char command[500]="";
		
	char buffer[BUFSIZE], buffer_2[BUFSIZE];

    char fdf_files[10][BUFSIZE];

    char merged_fdf_file[300000], merged_fdf_file_c[300000];

	char *dummy_buffer0;
	
	char dummyS[230]="";

	char *temp, *temp2;
	char *temp_line[10000], *temp_line2;


	int i=0, j=0, k; 


	int xL[100], yL[100];

	double t;

	double x_scale=1, y_scale=1, z_scale=1;
	
	double Norm_dx, Norm_dy, Norm_dz;

	double diffx_mx, diffy_mx, diffz_mx;

	int N_bands, N_spin, N_bands_k ;
	int N_DOS, N_pDOS;


	int ON_E_pDOS=0;
	int i_E_pDOS=0;

	int pDOS_n, pDOS_l,pDOS_m,pDOS_z;

	double pDOS_px, pDOS_py,pDOS_pz;

	char pDOS_spe[5]="";

	int ai, old_ai;

	int save_ON;

	int ON_newfile=0;

	int ON_read_Lattice=0;

	int orbit_i;

	int F_Band=0;

	int dummy_index=0, dummy_index2=0, dummy_index3=0;

	int N_read_Lattice=0;

	int N_LL=0;



	char str0[50]="", str1[50]="", str2[50]="", str3[50]="", str4[50]="";

	double L_x, L_y, L_z;

	double E_F;

    int N_fdf_file ;

	int N_i;

	int N_file =0;


////---- read BONDS_FINAL --------------------



  	//printf("args[1]= %s\n", args[1]);

    sprintf(command,"find %s -maxdepth 1 -name '*.BONDS_FINAL'", args[1]); 

    read_data0 = popen(command, "r");

    if(NULL == read_data0)  {  perror("There is no BONDS_FINAL");  return -1; }
	
    while (fgets(buffer, sizeof(buffer), read_data0))    {    N_file ++; }

	pclose(read_data0);

	if(N_file==0) {
		printf("There is no BONDS_FINAL\n");
		return -1;
	}


////---- read fdf files --------------------



  	sprintf(command,"find %s -maxdepth 1 -name '*.fdf'", args[1]); 

    read_data0 = popen(command, "r");

    if(NULL == read_data0)  {  perror("popen 실패");  return -1; }

    N_fdf_file =0 ;

    while (fgets(buffer, sizeof(buffer), read_data0))
    {

       
        Eliminate(buffer, '\r'); Eliminate(buffer, '\n'); 

        sprintf(fdf_files[N_fdf_file],  "%s", buffer);  	

        N_fdf_file ++;
    }

    pclose(read_data0);



	if (N_fdf_file ==0) 
	{
		printf("There is no fdf file.");
		return 0;
	} 

	
    for(i=0; i< N_fdf_file; i++ )
    {

        read_data0 = fopen(fdf_files[i], "rt");    

       	while( 1) 
    	{				
	        fgets( buffer, sizeof(buffer), read_data0 );

			if(feof( read_data0 )) break; 	
			
            sprintf(merged_fdf_file  ,  "%s%s", merged_fdf_file  , buffer);  	
			sprintf(merged_fdf_file_c,  "%s%s", merged_fdf_file_c, buffer);  	

        }    

        fclose(read_data0);    // 파일 포인터 닫기

    }



////----fdf �κ��� SystemLabel ����--------------------


    temp_line[0] = strtok(merged_fdf_file, "\n");     

	
	N_i=0; 

    while (temp_line[N_i] != NULL)               
    {
        
		N_i++;

        temp_line[N_i] = strtok(NULL, "\n");
	
		if(temp_line[N_i] == NULL) { break; }

		Eliminate(temp_line[N_i], '\r'); Eliminate(temp_line[N_i], '\n');			
			
    	temp = replaceAll(temp_line[N_i], "\t", " ")  ;  sprintf(temp_line[N_i],"%s",temp);		
		temp = replaceAll(temp_line[N_i], "	" , " ")  ;  sprintf(temp_line[N_i],"%s",temp);
			
    }


	
i=0;

for ( i =0 ; i<N_i; i++)
{	



		if(ON_read_Lattice==0)
		{   temp2=strtok(temp_line[i]," "); sprintf(str1, "%s", temp2); 
			if(strcmp(str1, "SystemLabel"               )==0) { temp2 = strtok(NULL  ," ");   sprintf(SystemLabel, "%s", temp2); }
			if(strcmp(str1, "NumberOfAtoms"             )==0) { temp2 = strtok(NULL  ," ");   N_atoms = atoi(temp2); }
			if(strcmp(str1, "SaveRho"                   )==0) { temp2 = strtok(NULL  ," ");   sprintf(On_save_R  , "%s", temp2); }
			if(strcmp(str1, "SaveElectrostaticPotential")==0) { temp2 = strtok(NULL  ," ");   sprintf(On_save_P  , "%s", temp2); }
			if(strcmp(str1, "%block"                    )==0) { temp2 = strtok(NULL  ," ");   
																								if(strcmp(temp2, "LatticeVectors"            )==0) {ON_read_Lattice=1;} 
																								if(strcmp(temp2, "BandLines"                 )==0) {F_Band         =1;}
																								if(strcmp(temp2, "PDOS.kgrid_Monkhorst_Pack" )==0) {F_DOS++          ;}
																								if(strcmp(temp2, "ProjectedDensityOfStates"  )==0) {F_DOS++          ;}
        														}
		}
		else      
		{ 
			temp2 = strtok(temp_line[i]," "); L_x=absol(atof(temp2)) ; if(L_x>0.00000001) { N_LL++;}
			temp2 = strtok(NULL  ," ")      ; L_y=absol(atof(temp2)) ; if(L_y>0.00000001) { N_LL++;}
			temp2 = strtok(NULL  ," ")      ; L_z=absol(atof(temp2)) ; if(L_z>0.00000001) { N_LL++;}
			N_read_Lattice++                ; if(N_read_Lattice==3) ON_read_Lattice=0; 
			
		}


}



if(N_LL!=3) { sprintf(On_save_R  , "F"); sprintf(On_save_P  , "F");}
	
//////////------------ ---- ������ġ ���� ���� ------------------------------------------------

Save_Bonds(args[1]);


VR_info(args[1]);


printf("where?--1-\n")	;


Save_S_V_R(args[1], merged_fdf_file_c);

printf("where?--2-\n")	;

if(strcmp(On_save_P, "T")==0 || strcmp(On_save_P, ".true.")==0) { Save_V(args[1]);   }
if(strcmp(On_save_R, "T")==0 || strcmp(On_save_R, ".true.")==0) { Save_Rho(args[1]); }
	
printf("where?--3-\n")	;

////// 
////////	/////------------ ---- Band data -> oneD ------------------------------------------------


if(F_Band==1)
{
    	sprintf(command,"/home3/DB_devel/Utils/new.gnubands < %s/%s.bands > %s/%s_out.band", args[1], SystemLabel, args[1], SystemLabel); system(command);	
		sprintf(filename2, "%s/%s_out.band", args[1], SystemLabel);  
		
		read_data2 = fopen(filename2,"rt");  

		while( 1 )
		{
				fgets( buffer, sizeof(buffer), read_data2 );

				if(buffer[0] !='#') { break;	}

				if(strcmp(buffer, "#\n")==0) { continue;	}
				
				temp=strtok(buffer," ");  
				temp=strtok(NULL  ," ");  

				if(strcmp(temp, "Nbands,")==0) 
				{
					temp=strtok(NULL  ," ");  
					temp=strtok(NULL  ," ");  
					temp=strtok(NULL  ," ");  
					temp=strtok(NULL  ," ");  N_bands=atoi(temp); 
					temp=strtok(NULL  ," ");  N_spin=atoi(temp); 
					temp=strtok(NULL  ," ");  N_bands_k=atoi(temp); 

				}

		}

		i=0; j=0;

		
		while( 1 )
		{
				fgets( buffer, sizeof(buffer), read_data2 );
				if(feof( read_data2 )) break;
		

				if(strcmp(buffer, "\n")==0) 
				{ 
					if(j!=0) i++; 
					j=0;  
					continue;
				}


				//Band_kE
				temp=strtok(buffer," ");  Band_kE[i][j][0] = atof(temp) ;
				temp=strtok(NULL  ," ");  Band_kE[i][j][1] = atof(temp) ;
				temp=strtok(NULL  ," ");  Band_spin[i][j]  = atoi(temp) ;

				j++;			
				
		}
		
	fclose(read_data2);

	
	
////	//////--------oneD plot writing -----------------------------
	
	sprintf(filename,  "%s/Bands.oneD", args[1]);  

	  fw = fopen(filename, "wt");
		
	////////------------------------------------------------------------------------------------

	   sprintf(dummyS,"#NumField: %d\n",     N_bands );	fprintf(fw, dummyS );

		   sprintf(dummyS,"#LabelX: k, LabelY: E (eV) \n");	fprintf(fw, dummyS );

		for(i=0; i<N_bands ; i++)
		{

		   sprintf(dummyS,"#Field%d: subband%d, NumPoints: %d \n",i,i, N_bands_k-1 );	fprintf(fw, dummyS );
		   for(j=0; j<N_bands_k-1 ; j++)   {	   fprintf(fw, "%lf, %lf\n", Band_kE[i][j][0], Band_kE[i][j][1]-E_F );	    }
		}

		  fclose(fw);       

		  // E_F shift �����ϵ�---------------------------------------

	sprintf(filename,  "%s/%s.bands", args[1], SystemLabel);  
	read_data0 = fopen(filename,"rt");  
	fgets( buffer, sizeof(buffer), read_data0 );		
		E_F=atof(buffer);
	fclose(read_data0);

	
	sprintf(filename,  "%s/info.dat", args[1]);  
	fw = fopen(filename,"wt");  

	fprintf(fw, "%s\n", SystemLabel );
	fprintf(fw, "%d\n", N_atoms );
	fprintf(fw, "%f\n", E_F );	
			
	fclose(fw);



}

////////////  /////------------ ---- DOS data -> oneD ------------------------------------------------

if(F_DOS==2)
{

    sprintf(filename2,  "%s/%s.DOS", args[1], SystemLabel);  
  
	read_data2 = fopen(filename2,"rt");  

	if( read_data2 != NULL)
	{
		
		i=0; 

		while( 1 )
		{
				fgets( buffer, sizeof(buffer), read_data2 );
				if(feof( read_data2 )) break;

				temp=strtok(buffer," ");  DOS[i][0] = atof(temp) ; 
				temp=strtok(NULL  ," ");  DOS[i][1] = atof(temp) ; 			
				i++;			
		}
		N_DOS = i;
		fclose(read_data2);
////	//--------oneD plot writing -----------------------------

	sprintf(filename,  "%s/DOS.oneD", args[1]);  

	  fw = fopen(filename, "wt");
		
	   sprintf(dummyS,"#NumField: %d\n",     1 )              ;	fprintf(fw, dummyS );
	   sprintf(dummyS,"#LabelX: k, LabelY: E (eV) \n")        ;	fprintf(fw, dummyS );
	   sprintf(dummyS,"#Field: DOS, NumPoints: %d \n", N_DOS );	fprintf(fw, dummyS );

	   i=0;


	   for(i=0; i<N_DOS ; i++)
	   {

		   fprintf(fw, "%lf, %lf\n", DOS[i][0], DOS[i][1] );	   
	   }

	  fclose(fw);    

////	/////------------ ---- pDOS data -> oneD ------------------------------------------------
////
////	  //--- read PDOS ----------------
////		//sprintf(filename2,  "%s", args[7]);  	
////
		sprintf(filename2,  "%s/%s.PDOS", args[1], SystemLabel);  
		
		read_data2 = fopen(filename2,"rt");  	
		while( 1 )
		{
				fgets( buffer, sizeof(buffer), read_data2 ); temp=strtok(buffer," ");   
				if(strcmp(temp, "<energy_values"    )==0) { i=0; }
				if(strcmp(temp, "</energy_values>\n")==0) { break;}
				i++;
		}
		fclose(read_data2);

		N_pDOS = i-1;

//	//--- �⺻���� �б�  ������ �� �а� ���� ----------------
		read_data2 = fopen(filename2,"rt");  	
		while( 1 )
		{

				fgets( buffer, sizeof(buffer), read_data2 ); temp=strtok(buffer," ");   
				if(feof( read_data2 )) break;
				
				if(ON_E_pDOS==1) 
				{
					E_pDOS[i_E_pDOS]= atof(temp)-E_F;

					i_E_pDOS++;

					if(i_E_pDOS==N_pDOS) {ON_E_pDOS =0;}	
				}

				if(strcmp(temp, "<energy_values") ==0) ON_E_pDOS=1; 			
				
				if(strcmp(temp, "atom_index=\"")  ==0) 
				{				
					temp=strtok(NULL  ," "); Eliminate(temp, ' '); Eliminate(temp, '"'); N_orbit_per_atom[atoi(temp)] ++;
					
				}

		}
		fclose(read_data2);
	
		read_data2 = fopen(filename2,"rt");  	

		ai=0;
		old_ai=ai;
		orbit_i=0;
		i_E_pDOS=0;

		
		while( 1 )
		{
			fgets( buffer, sizeof(buffer), read_data2 ); temp=strtok(buffer," ");   
			if(feof( read_data2 )) break;

			if(strcmp(temp, "</data>\n")      ==0) 
			{ 
				save_ON = 0;  
				if(orbit_i==N_orbit_per_atom[ai]) {	
					//fclose(fw);

					sprintf(dummyS,"%s/pDOS_atom%d_%s_x=%.2lf_y=%.2lf_z=%.2lf.oneD", args[1], ai, pDOS_spe, pDOS_px, pDOS_py, pDOS_pz  );
					fwT = fopen(dummyS, "wt");				
					fprintf(fwT, "#NumField: 1\n" );
					fprintf(fwT, "#LabelX: E (eV), LabelY: pDOS \n");
					fprintf(fwT, "#Field%d: Total/atom, NumPoints: %d \n", 1, N_pDOS);

				    for(int p_i=0; p_i<N_pDOS; p_i++)
					{
						fprintf(fwT, "%lf, %lf\n", pDOS_T[p_i], pDOS_E_T[p_i]);
					}
				
					fclose(fwT);

				    for(int p_i=0; p_i<N_pDOS; p_i++)
					{
						pDOS_E_T[p_i] = 0.0;
					}



				}
				
			}

			if(save_ON==1) 
			{

				//fprintf(fw, "%lf, %lf\n", E_pDOS[i_E_pDOS++], atof(temp));			

				pDOS_T[i_E_pDOS]   = E_pDOS[i_E_pDOS];
				pDOS_E_T[i_E_pDOS] += atof(temp);

				i_E_pDOS++;

			}

			if(strcmp(temp, "<data>\n")       ==0) 
			{
				if(ON_newfile==1)
				{
					//sprintf(dummyS,"%s/pDOS_atom%d_%s_x=%lf_y=%lf_z=%lf.oneD", args[1], ai, pDOS_spe, pDOS_px, pDOS_py, pDOS_pz  );
					//fw = fopen(dummyS, "wt");				
					//fprintf(fw, "#NumField: %d\n", N_orbit_per_atom[ai] );
					//fprintf(fw, "#LabelX: E (eV), LabelY: pDOS \n");

					//sprintf(dummyS,"%s/pDOS_atom%d_%s_x=%.2lf_y=%.2lf_z=%.2lf.oneD", args[1], ai, pDOS_spe, pDOS_px, pDOS_py, pDOS_pz  );
					//fwT = fopen(dummyS, "wt");				
					//fprintf(fwT, "#NumField: 1\n" );
					//fprintf(fwT, "#LabelX: E (eV), LabelY: pDOS \n");
					//fprintf(fwT, "#Field%d: Total/atom, NumPoints: %d \n", 1, N_pDOS);

					ON_newfile=0;

				}

				orbit_i++; 
				//fprintf(fw, "#Field%d: (nlmz)=(%d%d%d%d), NumPoints: %d \n", orbit_i, pDOS_n, pDOS_l, pDOS_m, pDOS_z, N_pDOS);

				save_ON =1 ;
				i_E_pDOS=0 ;
			}

			if(strcmp(temp, "atom_index=\"")  ==0) 
			{ 
				old_ai = ai;

				temp=strtok(NULL  ," "); Eliminate(temp, ' '); Eliminate(temp, '"'); ai     = atoi(temp); 

				if(ai!=old_ai)	{ ON_newfile=1;		orbit_i=0;	 }
			}
			if(strcmp(temp, "n=\"")           ==0) { temp=strtok(NULL  ," "); Eliminate(temp, ' '); Eliminate(temp, '"'); pDOS_n = atoi(temp); 	}
			if(strcmp(temp, "l=\"")           ==0) { temp=strtok(NULL  ," "); Eliminate(temp, ' '); Eliminate(temp, '"'); pDOS_l = atoi(temp); 	}
			if(strcmp(temp, "m=\"")           ==0) { temp=strtok(NULL  ," "); Eliminate(temp, ' '); Eliminate(temp, '"'); pDOS_m = atoi(temp); 	}
			if(strcmp(temp, "z=\"")           ==0) { temp=strtok(NULL  ," "); Eliminate(temp, ' '); Eliminate(temp, '"'); pDOS_z = atoi(temp); 	}

			if(strcmp(temp, "position=\"")    ==0) { temp=strtok(NULL  ," "); Eliminate(temp, '"'); pDOS_px= atof(temp); 
			                                         temp=strtok(NULL  ," "); Eliminate(temp, '"'); pDOS_py= atof(temp); 
													 temp=strtok(NULL  ," "); Eliminate(temp, '"'); pDOS_pz= atof(temp); 
			                                       }

            if(strstr(temp, "species")     !=NULL) { temp=strtok(temp,"\"") ; temp=strtok(NULL,"\"");	sprintf(pDOS_spe,temp);	}
		}

		fclose(read_data2);

	}

}




}
// ���� ��ȣ ��ȣ ����--------
int Num_atom(char*  atom_string)
{

	char   Name_Atom[119][3]={"0" ,
		                      "H" ,"He","Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar",
		                      "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
							  "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
							  "Cs","Ba",
							  "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
							  "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
							  "Fr","Ra",
							  "Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
							  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
	
	};	

	int i;

   char *ptr;

   char *temp;
   char temp2[10];
   char *atom_st; 
   
   strcpy(temp2,atom_string);
 
   ptr = index( temp2, '-');

	if ( ptr )   { temp=strtok(temp2,"-");   }
	else{ temp=temp2 ; }

	 
for(i=1;i<=118; i++) if(strcmp(temp, Name_Atom[i])==0) { return i;}

	return 0;
		
}
//------------------------------------------------------
void Save_Bonds(char* path)
{
	FILE *read_data0, *read_data2;
	char filename[230]="";	
	char command[500]="";
	char *temp, *temp2;
	char buffer[BUFSIZE], buffer_2[BUFSIZE];

	int i, j;

	char temp_atom[20]="";	

	

	sprintf(filename,  "%s/%s.BONDS_FINAL", path, SystemLabel);  

	read_data0 = fopen(filename,"rt");

	i=0;

	while( 1 )
	{
		
	        fgets( buffer, sizeof(buffer), read_data0 );	
			
			if(feof( read_data0 )) break;

			sprintf(buffer, "%s", replaceAll(buffer, "-", " -"));				
			
			temp=strtok(buffer," ");  
			
			if(strcmp(temp, "Neighbors")==0) 
			{
			
				temp = strtok(NULL  ," "); 
				temp = strtok(NULL  ," ");  sprintf(ID_Atom[i], "%s" , temp); sprintf(ID_cAtom[c_i], "%s", temp); 
				temp = strtok(NULL  ," ");  sprintf(temp_atom,"%s",temp); 
				temp = strtok(NULL  ," ");  
				temp = strtok(NULL  ," ");  P_Atom[i][1] = atof(temp); 
				temp = strtok(NULL  ," ");  
				if ( point_count(temp) == 1) 
				{ 
					P_Atom[i][2] = atof(temp);
					temp = strtok(NULL  ," ");  P_Atom[i][3] = atof(temp);
				} 
				else if( point_count(temp) == 2)
				{
					char dummy[50]="00";
					char * point_p ;

					sprintf(dummy,"00%s",temp);
					point_p = strrchr(temp,'.')-3;
					P_Atom[i][3] = atof(point_p);
					*point_p = '\0';
					P_Atom[i][2] = atof(temp);
				}

				N_P_Atom[i]=Num_atom(temp_atom);  
				sprintf(ID_Atom[i],  "%s%s %d", ID_Atom[i], temp_atom, N_P_Atom[i] ); 


				c_i++;
				cL_i=i;

			

			}			
			else
			{				
				
				                            sprintf(ID_Atom[i],  "%s" , temp);  sprintf(Link[i], "%d %d", cL_i, i);  
				
				temp = strtok(NULL  ," ");  sprintf(temp_atom,"%s",temp);  
				temp = strtok(NULL  ," "); 
				temp = strtok(NULL  ," "); 
				temp = strtok(NULL  ," "); 
				temp = strtok(NULL  ," "); 
				temp = strtok(NULL  ," "); P_Atom[i][1] = atof(temp);  
				temp = strtok(NULL  ," "); 

				if ( point_count(temp) == 1) 
				{ 
					P_Atom[i][2] = atof(temp);
					temp = strtok(NULL  ," ");  P_Atom[i][3] = atof(temp);
				} 
				else if( point_count(temp) == 2)
				{
					char dummy[50]="00";
					char * point_p ;

					sprintf(dummy,"00%s",temp);

					point_p = strrchr(temp,'.')-3;
					P_Atom[i][3] = atof(point_p);

					*point_p = '\0';

					P_Atom[i][2] = atof(temp);
				}
				

				N_P_Atom[i]=Num_atom(temp_atom);  


				sprintf(ID_Atom[i],  "%s%s %d", ID_Atom[i], temp_atom, N_P_Atom[i]); 


	
			}

			i++;

	}

	fclose(read_data0);

	N_V_atoms = i;






	for(i=0; i<N_V_atoms; i++)
	{

		Is_same=0;	

		for(j=i; j<N_V_atoms; j++)
		{
			if(i==j) continue;
		
			if(P_Atom[i][1]==P_Atom[j][1] && P_Atom[i][2]==P_Atom[j][2] && P_Atom[i][3]==P_Atom[j][3]) { Is_same=1; OverLap_Atom[i]=1; break; }
		}
		
	}


	
	N_Links = N_V_atoms-c_i;


}
//------------------------------------------------------
void VR_info(char* path)
{
	
	FILE *read_data0, *read_data2;
	FILE *fr, *fw, *fr_data, *fw_data;
	FILE *fw_SVR;

	char command[500]="";
	char *temp, *temp2;
	char buffer[BUFSIZE], buffer_2[BUFSIZE];

	char filename[230]="";	



	sprintf(filename,"%s/%s.VH",path,SystemLabel);
	if (read_data0 = fopen(filename, "r")) { fclose(read_data0);   }
    else                                   { sprintf(On_save_P  , "F"); printf("No VH file-----------------------------\n");		    }



	if(strcmp(On_save_P, "T")==0 || strcmp(On_save_P, ".true.")==0) 
	{ 		
		
		sprintf(filename,"%s/input_vh.inp",path);
		

		fw = fopen(filename, "wt"); 

		fprintf(fw, SystemLabel );	

		fprintf(fw, "\nvh\n" );	

		fprintf(fw, "0.0 0.0 0.0\n" );	

		fprintf(fw, "1\n" );	

		fprintf(fw, "unformatted" );	 

		fclose(fw);



		sprintf(command,"cd %s;/home3/SWs/SIESTA/siesta-master/Util/Grid/g2c_ng -x 1 -y 1 -z 1 -n 1 -s %s.STRUCT_OUT -g %s.VH", path, SystemLabel, SystemLabel); system(command);	


		sprintf(command,"mv %s/Grid.cube %s/%s.VH.cube", path, path, SystemLabel); system(command);

		sprintf(filename,  "%s/%s.VH.cube", path, SystemLabel);  	

		read_data0 = fopen(filename,"rt");  

		fgets( buffer, sizeof(buffer), read_data0 );	
		fgets( buffer, sizeof(buffer), read_data0 );	
		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  N_atoms = atoi(temp) ; 			

		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  Nx = atoi(temp) ; temp=strtok(NULL  ," "); dx = atof(temp) ; temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); 
		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  Ny = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," "); dy = atof(temp); temp=strtok(NULL  ," "); 
		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  Nz = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); dz = atof(temp) ; 

		Nzy = Nz*Ny;

		N_nod = Nx*Ny*Nz ;

		fclose(read_data0);


	}



	

	sprintf(filename,"%s/%s.RHO",path,SystemLabel);
	if (read_data0 = fopen(filename, "r")) { fclose(read_data0);   }
    else                                   { sprintf(On_save_R  , "F"); printf("No RHO file-----------------------------\n");		    }

	if(strcmp(On_save_R, "T")==0 || strcmp(On_save_R, ".true.")==0) 
	{ 

	
		sprintf(filename,"%s/input_rho.inp", path);
		fw = fopen(filename, "wt"); fprintf(fw, SystemLabel );	fprintf(fw, "\nrho\n" );	fprintf(fw, "0.0 0.0 0.0\n" );	fprintf(fw, "1\n" );	fprintf(fw, "unformatted" );	 fclose(fw);
		
		sprintf(command,"cd %s;/home3/SWs/SIESTA/siesta-master/Util/Grid/g2c_ng -x 1 -y 1 -z 1 -n 1 -s %s.STRUCT_OUT -g %s.RHO", path, SystemLabel, SystemLabel); system(command);	
		sprintf(command,"mv %s/Grid.cube %s/%s.RHO.cube", path, path, SystemLabel); system(command);

		sprintf(filename,  "%s/%s.RHO.cube", path, SystemLabel);  	

		read_data0 = fopen(filename,"rt");  

		fgets( buffer, sizeof(buffer), read_data0 );	
		fgets( buffer, sizeof(buffer), read_data0 );	
		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  N_atoms = atoi(temp) ; 			

		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  Nx = atoi(temp) ; temp=strtok(NULL  ," "); dx = atof(temp) ; temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); 
		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  Ny = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," "); dy = atof(temp); temp=strtok(NULL  ," "); 
		fgets( buffer, sizeof(buffer), read_data0 ); temp=strtok(buffer," ");  Nz = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); dz = atof(temp) ; 

		Nzy = Nz*Ny;

		N_nod = Nx*Ny*Nz ;

		fclose(read_data0);

	}


	


}
//---------------------------------------
int Save_V(char* path)
{
	

	FILE *read_data0, *read_data2;
	FILE *fr, *fw, *fr_data, *fw_data;

	FILE *fw_SVR;

		char command[500]="";
char *temp, *temp2;
char buffer[BUFSIZE], buffer_2[BUFSIZE];

char filename[230]="";	

int i,j,k, m, order;

double xx, yy, zz;
char dummyS[230]="";

	xmx=-1e+100, ymx=-1e+100, zmx=-1e+100;
    xmn= 1e+100, ymn= 1e+100, zmn= 1e+100;
    

	//sprintf(filename,"input_vh.inp");

	//fw = fopen(filename, "wt"); fprintf(fw, SystemLabel );	fprintf(fw, "\nvh\n" );	fprintf(fw, "0.0 0.0 0.0\n" );	fprintf(fw, "1\n" );	fprintf(fw, "unformatted" );	 fclose(fw);

	////sprintf(command,"/SYSTEM/Siesta/4.0.1/Util/Grid/grid2cube < %s", filename); system(command);
	//sprintf(command,"/home3/SWs/SIESTA/siesta-master/Util/Grid/g2c_ng -x 1 -y 1 -z 1 -n 1 -s %s.STRUCT_OUT -g %s.VH", SystemLabel, SystemLabel); system(command);	

	//sprintf(command,"mv Grid.cube %s.VH.cube", SystemLabel); system(command);

	

	sprintf(filename,  "%s/%s.VH.cube", path, SystemLabel);  	

	//sprintf(filename,  "%s", args[3]);  
	
	read_data0 = fopen(filename,"rt");  

	fgets( buffer, sizeof(buffer), read_data0 );	
	fgets( buffer, sizeof(buffer), read_data0 );	
	fgets( buffer, sizeof(buffer), read_data0 ); //temp=strtok(buffer," ");  N_atoms = atoi(temp) ; 			

	fgets( buffer, sizeof(buffer), read_data0 ); //temp=strtok(buffer," ");  Nx = atoi(temp) ; temp=strtok(NULL  ," "); dx = atof(temp) ; temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); 
	fgets( buffer, sizeof(buffer), read_data0 ); //temp=strtok(buffer," ");  Ny = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," "); dy = atof(temp); temp=strtok(NULL  ," "); 
	fgets( buffer, sizeof(buffer), read_data0 ); //temp=strtok(buffer," ");  Nz = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); dz = atof(temp) ; 

	Nzy = Nz*Ny;

	N_nod = Nx*Ny*Nz ;

	if(N_nod > MAX_RhoV) 
	{
		sprintf(filename, "%s/Structure_V_R.js3D", path);

		fw_SVR = fopen(filename, "at");

		sprintf(dummyS,"@EP_-1 0\n" );

		fprintf(fw_SVR, dummyS);

		fclose(fw_SVR);    

		return 1;
	}

	//printf("Nx, Ny, Nz= %d %d %d \n",Nx, Ny, Nz);

	

	for(i=0 ; i<N_atoms ; i++)
	{
		fgets( buffer, sizeof(buffer), read_data0 );	
	}
	


	i=0;


	j=0;


	while( 1 )
	{
			fgets( buffer, sizeof(buffer), read_data0 );

			j++;

			
			if(feof( read_data0 )) {  break;}

			temp=strtok(buffer," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz; 
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i] > func_mx) func_mx = func[i]  ; if ( func[i]  < func_mn) func_mn = func[i];			   
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz; 			       
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func_mx) func_mx = func[i]  ; if ( func[i]  < func_mn) func_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz; 
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func_mx) func_mx = func[i]  ; if ( func[i]  < func_mn) func_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz;
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func_mx) func_mx = func[i]  ; if ( func[i]  < func_mn) func_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz;
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func_mx) func_mx = func[i]  ; if ( func[i]  < func_mn) func_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz;
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func_mx) func_mx = func[i]  ; if ( func[i]  < func_mn) func_mn = func[i];
			   i++;
			
	
	}

	
	//printf("Nx, Ny, Nz= %f %f %f \n",xmx,ymx);

		

fclose(read_data0);

sprintf(filename, "%s/Structure_V_R.js3D", path);


fw_SVR = fopen(filename, "at");



		sprintf(dummyS,"@EP_data %d\n", N_nod );

		fprintf(fw_SVR, dummyS);

		for (ti=0;ti<N_t;ti++)
		{		

			for(i=0; i<Nx; i++) for(j=0; j<Ny; j++) 
			{
				
				for(k=0; k<Nz; k++)
				{
					double Norm_func;
					order = k + j*Nz + i*Nzy ;			
					Norm_func=(func[order]-func_mn)/(func_mx-func_mn);	    	
				   fprintf(fw_SVR, "%.2e\n", Norm_func );

				}
			}
		}


    	fprintf(fw_SVR, "\n");
	


fclose(fw_SVR);      

return 0;



}
//------------------------------------------------------
int Save_Rho(char* path)
{
		FILE *read_data0, *read_data2;
	FILE *fr, *fw, *fr_data, *fw_data;
	FILE *fw_SVR;

		char command[500]="";
char *temp, *temp2;
char buffer[BUFSIZE], buffer_2[BUFSIZE];

char filename[230]="", filename2[230]="";	

int i,j,k, m, order;

double xx, yy, zz;

char dummyS[230]="";

	xmx=-1e+100, ymx=-1e+100, zmx=-1e+100;
    xmn= 1e+100, ymn= 1e+100, zmn= 1e+100;



	/*sprintf(filename,"input_rho.inp");

	fw = fopen(filename, "wt"); fprintf(fw, SystemLabel );	fprintf(fw, "\nrho\n" );	fprintf(fw, "0.0 0.0 0.0\n" );	fprintf(fw, "1\n" );	fprintf(fw, "unformatted" );	 fclose(fw);
	
	sprintf(command,"/home3/SWs/SIESTA/siesta-master/Util/Grid/g2c_ng -x 1 -y 1 -z 1 -n 1 -s %s.STRUCT_OUT -g %s.RHO", SystemLabel, SystemLabel); system(command);	

	sprintf(command,"mv Grid.cube %s.RHO.cube", SystemLabel); system(command);*/

	sprintf(filename2, "%s/%s.RHO.cube", path, SystemLabel);





		//sprintf(filename2,  "%s.RHO.cube", SystemLabel);  
		
		read_data2 = fopen(filename2,"rt");  


		fgets( buffer, sizeof(buffer), read_data2 );  
		fgets( buffer, sizeof(buffer), read_data2 ); 
		fgets( buffer, sizeof(buffer), read_data2 ); //temp=strtok(buffer," ");  N_atoms = atoi(temp) ; 	
		fgets( buffer, sizeof(buffer), read_data2 ); //temp=strtok(buffer," ");  Nx = atoi(temp) ; temp=strtok(NULL  ," "); dx = atof(temp) ; temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); 
		fgets( buffer, sizeof(buffer), read_data2 ); //temp=strtok(buffer," ");  Ny = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," "); dy = atof(temp); temp=strtok(NULL  ," "); 
		fgets( buffer, sizeof(buffer), read_data2 ); //temp=strtok(buffer," ");  Nz = atoi(temp) ; temp=strtok(NULL  ," ");                   temp=strtok(NULL  ," ");                  temp=strtok(NULL  ," "); dz = atof(temp) ; 

		Nzy = Nz*Ny;



		N_nod = Nx*Ny*Nz ;


		if(N_nod > MAX_RhoV) 
		{
			sprintf(filename, "%s/Structure_V_R.js3D", path);

				

			fw_SVR = fopen(filename, "at");
			sprintf(dummyS,"@Rho_-1 0\n" );
			fprintf(fw_SVR, dummyS);
			fclose(fw_SVR);    
			return 1;
		}



		for(i=0 ; i<N_atoms ; i++)
		{
			fgets( buffer, sizeof(buffer), read_data2 );	
		}

		i=0;

		j=0;
		
		while( 1 )
		{
			j++;
				fgets( buffer, sizeof(buffer), read_data2 );

				if(feof( read_data2 )) break;


			temp=strtok(buffer," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz; 
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i] > func2_mx) func2_mx = func[i]  ; if ( func[i]  < func2_mn) func2_mn = func[i];			   
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz; 			       
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func2_mx) func2_mx = func[i]  ; if ( func[i]  < func2_mn) func2_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz; 
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func2_mx) func2_mx = func[i]  ; if ( func[i]  < func2_mn) func2_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz;
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func2_mx) func2_mx = func[i]  ; if ( func[i]  < func2_mn) func2_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz;
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func2_mx) func2_mx = func[i]  ; if ( func[i]  < func2_mn) func2_mn = func[i];
			   i++;

			temp=strtok(NULL  ," ");  if(temp==NULL) { continue; } func[i] = atof(temp) ;   xi =i/Nzy; md = i%Nzy ; yi = md/Nz; zi = md%Nz; xx=xi*dx; yy=yi*dy; zz=zi*dz;
	           if ( xx > xmx ) xmx = xx; if ( xx < xmn ) xmn = xx; if ( yy > ymx ) ymx = yy; if ( yy < ymn ) ymn = yy; if ( zz > zmx ) zmx = zz; if ( zz < zmn ) zmn = zz; if ( func[i]   > func2_mx) func2_mx = func[i]  ; if ( func[i]  < func2_mn) func2_mn = func[i];
			   i++;

				
		}
		
		
	fclose(read_data2);

		


	sprintf(filename, "%s/Structure_V_R.js3D", path);

	fw_SVR = fopen(filename, "at");

		

	sprintf(dummyS,"@Rho_data %d\n", N_nod );

	fprintf(fw_SVR, dummyS);

	for (ti=0;ti<N_t;ti++)
	{		
	

		for(i=0; i<Nx; i++) for(j=0; j<Ny; j++) 
		{
			for(k=0; k<Nz; k++)
			{

			double Norm_func;

			order = k + j*Nz + i*Nzy ;			

			Norm_func=(func[order]-func2_mn)/(func2_mx-func2_mn);			
	   	   
		   fprintf(fw_SVR, "%.2e\n", Norm_func );					
				   

				   
			}
		}


	}


	fprintf(fw_SVR, "\n");

				
	

	//
	fclose(fw_SVR);       

	return 0;



}

//------------------------------------------------------
void Save_S_V_R(char* path, char* fdf_file_contents)
{

	FILE *read_data0, *read_data2;
	FILE *fr, *fw, *fr_data, *fw_data;

	FILE *fw_SVR;

	char command[500]="";
	char *temp, *temp2;
	char buffer[BUFSIZE], buffer_2[BUFSIZE];
	char filename[2300]="", filename2[2300]="";	
	int i,j,k, m, order;

	int tti;

	int fi;

	//int dum_ON=0;

	char dummyS[10000]="";
	char dummyS2[10000]="";
	

	i=0;

	


    sprintf(filename, "%s/Structure_V_R.js3D", path);
	fw_SVR = fopen(filename, "wt");

	sprintf(dummyS, "@params 18\n");	
	sprintf(dummyS, "%sNx=%d\n"    , dummyS, Nx);
	sprintf(dummyS, "%sNy=%d\n"    , dummyS, Ny);
	sprintf(dummyS, "%sNz=%d\n"    , dummyS, Nz);
	sprintf(dummyS, "%sN_nod=%d\n" , dummyS, N_nod);
	sprintf(dummyS, "%sN_t=%d\n"   , dummyS, N_t);
	sprintf(dummyS, "%sxmx=%lf\n"  , dummyS, xmx);
	sprintf(dummyS, "%symx=%lf\n"  , dummyS, ymx);
	sprintf(dummyS, "%szmx=%lf\n"  , dummyS, zmx);

	sprintf(dummyS, "%sxmn=%lf\n"  , dummyS, xmn);
	sprintf(dummyS, "%symn=%lf\n"  , dummyS, ymn);
	sprintf(dummyS, "%szmn=%lf\n"  , dummyS, zmn);

    sprintf(dummyS, "%sdx=%lf\n"  , dummyS, dx);
	sprintf(dummyS, "%sdy=%lf\n"  , dummyS, dy);
	sprintf(dummyS, "%sdz=%lf\n"  , dummyS, dz);
    
	sprintf(dummyS, "%sP_max=%lf\n"  , dummyS, func_mx);
	sprintf(dummyS, "%sP_min=%lf\n"  , dummyS, func_mn);	
	
	sprintf(dummyS, "%sR_max=%lf\n"  , dummyS, func2_mx);
	sprintf(dummyS, "%sR_min=%lf\n"  , dummyS, func2_mn);	

	temp = replaceAll(fdf_file_contents, "%block", "%%block")  ;  sprintf(fdf_file_contents,"%s",temp);		
	temp = replaceAll(fdf_file_contents, "%endblock", "%%endblock")  ;  sprintf(fdf_file_contents,"%s",temp);		
	
    strcat(dummyS, fdf_file_contents);        // s2 뒤에 s1를 붙임

	strcat(dummyS, "\n");
	
	fprintf(fw_SVR, dummyS );

	printf("fdf_file_contents\n %s\n", dummyS);
       
	
////////-------- ���� ��ġ ���� ----------------------------------------------------------------------------

    sprintf(dummyS,"@Atom_coor %d\n", N_V_atoms ); fprintf(fw_SVR, dummyS );	

	for(i=0; i<N_V_atoms; i++)
	{

////		sprintf(dummyS,"%s %lf %lf %lf %d\n", ID_Atom[i], P_Atom[i][1], P_Atom[i][2], P_Atom[i][3], OverLap_Atom[i] );	fprintf(fw_SVR, dummyS );
		sprintf(dummyS,"%s %d %lf %lf %lf %d\n", ID_Atom[i], N_P_Atom[i], P_Atom[i][1], P_Atom[i][2], P_Atom[i][3], OverLap_Atom[i] );	
		fprintf(fw_SVR, dummyS );	
	}

	    sprintf(dummyS,"@Bond_coor %d\n", N_Links );	
		fprintf(fw_SVR, dummyS );	
		
	for(i=0; i<N_V_atoms; i++)
	{		


		if(strcmp(Link[i], "")==0) { continue;}
		sprintf(dummyS,"%s %d\n", Link[i], OverLap_Atom[i] );	
		fprintf(fw_SVR, dummyS );	
	
	}



	//fprintf(fw_SVR, dummyS );

	fclose(fw_SVR);       

	

	//sprintf(dummyS2,"%s/Structure_V_R.js3D", path );	

	//fw = fopen(dummyS2, "wt");
		   
	//fprintf(fw, "%s", dummyS);
	
    //fclose(fw); 



	//fw = fopen("path.txt", "wt");
		   
	//fprintf(fw, "%s", path);
	
    //fclose(fw); 

}

int orderF(int i, int j, int k){	return k + j*Nz + i*Nzy ; }
//------------------------------------------------------

// ���ڿ� ����----------

void Eliminate(char *str, char ch)
{

    for (; *str != '\0'; str++)//���� ���ڸ� ���� ������ �ݺ�

    {

        if (*str == ch)//ch�� ���� ������ ��

        {

            strcpy(str, str + 1);

            str--;

        }

    }

}


//���ڿ� ġȯ-------------------------

char *replaceAll(char *s, const char *olds, const char *news) 
{
  char *result, *sr;
  size_t i, count = 0;
  size_t oldlen = strlen(olds); if (oldlen < 1) return s;
  size_t newlen = strlen(news);


  if (newlen != oldlen) {
    for (i = 0; s[i] != '\0';) {
      if (memcmp(&s[i], olds, oldlen) == 0) count++, i += oldlen;
      else i++;
    }
  } else i = strlen(s);


  result = (char *) malloc(i + 1 + count * (newlen - oldlen));
  if (result == NULL) return NULL;


  sr = result;
  while (*s) {
    if (memcmp(s, olds, oldlen) == 0) {
      memcpy(sr, news, newlen);
      sr += newlen;
      s  += oldlen;
    } else *sr++ = *s++;
  }
  *sr = '\0';

  return result;
}


//----------------------- ���밪 ���----------------------------------
double absol(double value)
{

    if(value>=0.0)
    {
        return value;

    }else
    {
        return -value;
    }



}

//// ���ڿ� ���� ���鹮�� ���� �Լ�
char* rtrim(char* s) {
  char t[MAX_STR_LEN];
  char *end;

  // Visual C 2003 ���Ͽ�����
  // strcpy(t, s);
  // �̷��� �ؾ� ��
  strcpy(t, s); // �̰��� Visual C 2005��
  end = t + strlen(t) - 1;
  while (end != t && isspace(*end))
    end--;
  *(end + 1) = '\0';
  s = t;

  return s;
}
//
//
//// ���ڿ� ���� ���鹮�� ���� �Լ�
char* ltrim(char *s) {
  char* begin;
  begin = s;

  while (*begin != '\0') {
    if (isspace(*begin))
      begin++;
    else {
      s = begin;
      break;
    }
  }

  return s;
}


// ���ڿ� �յ� ���� ��� ���� �Լ�
char* trim(char *s) {
  return rtrim(ltrim(s));
}

int point_count(char * parent_string )
{
    int count = 0;    
    int exampleLength = strlen(parent_string) ; // size of array
	//printf("exampleLength %s %d \n", parent_string, exampleLength);
    for (int i = 0; i < exampleLength; i++) {
        if ( parent_string[i]  == '.') count++;
    }
    return count;
}

