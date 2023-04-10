import time
import os
from pymongo import MongoClient

#------------------------------------------------------------------------------------------		
def atom_info_html_table():	  
	#pt_atom_num = [[0]*18]*10
	pt_atom_num = [[0 for _ in range(18)] for _ in range(10)]		 

	atom_p_pt = [ 
	[0,0],																										  [0,17],
	[1,0],[1,1],																[1,12],[1,13],[1,14],[1,15],[1,16], [1,17],
	[2,0],[2,1],																[2,12],[2,13],[2,14],[2,15],[2,16], [2,17],
	[3,0],[3,1],[3,2], [3,3],[3,4],[3,5],[3,6],[3,7],[3,8],[3,9],[3,10],[3,11], [3,12],[3,13],[3,14],[3,15],[3,16], [3,17],
	[4,0],[4,1],[4,2], [4,3],[4,4],[4,5],[4,6],[4,7],[4,8],[4,9],[4,10],[4,11], [4,12],[4,13],[4,14],[4,15],[4,16], [4,17],
	[5,0],[5,1],[5,2], 
					   [7,3],[7,4],[7,5],[7,6],[7,7],[7,8],[7,9],[7,10],[7,11],[7,12],[7,13],[7,14],[7,15],[7,16],
					   [5,3],[5,4],[5,5],[5,6],[5,7],[5,8],[5,9],[5,10],[5,11],[5,12],[5,13],[5,14],[5,15],[5,16],[5,17],
	[6,0],[6,1],[6,2], 
					   [8,3],[8,4],[8,5],[8,6],[8,7],[8,8],[8,9],[8,10],[8,11],[8,12],[8,13],[8,14],[8,15],[8,16],
					   [6,3],[6,4],[6,5],[6,6],[6,7],[6,8],[6,9],[6,10],[6,11],[6,12],[6,13],[6,14],[6,15],[6,16],[6,17]
	]
   
	c_i = 1
	for atom_p in atom_p_pt:
		pt_atom_num[atom_p[0]][atom_p[1]] = c_i		
		c_i += 1
	
	Atom_Name = ["",
			   "H" ,																				 "He",
			   "Li","Be",												   "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
			   "Na","Mg",												   "Al","Si","P" ,"S" ,"Cl","Ar",
			   "K" ,"Ca","Sc", "Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
			   "Rb","Sr","Y" , "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
			   "Cs","Ba","La", 
							   "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
							   "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
			   "Fr","Ra","Ac", 
							   "Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
							   "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"  ]

	Atom_BG_color = ["000000",
			   "388e3c",																																				"fbc02d",
			   "fa9d00","e64a19",																						  "5c3f37","c1185b","c1185b","c1185b","f37b00","fbc02d",
			   "fa9d00","e64a19",																						  "00786a","5c3f37","c1185b","c1185b","f37b00","fbc02d",
			   "fa9d00","e64a19","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","00786a","5c3f37","5c3f37","c1185b","f37b00","fbc02d",
			   "fa9d00","e64a19","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","00786a","00786a","5c3f37","5c3f37","f37b00","fbc02d",
			   "fa9d00","e64a19","2f3e9c", 
										  "2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c","2f3e9c",
										  "0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","0094a4","00786a","00786a","00786a","00786a","5c3f37","fbc02d",
			   "fa9d00","e64a19","7a1f91", 
										  "7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91","7a1f91",
										  "0094a4","0094a4","606060","606060","606060","606060","606060","606060","606060","606060","606060","606060","606060","606060","606060"  ]


	pseudo_type_VASP = [  ['LDA','PBE'], 
	['LDA','PBE'],																																																								 ['LDA','PBE'],
	['LDA','PBE'],['LDA','PBE'],																																			 ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
	['LDA','PBE'],['LDA','PBE'],																																			 ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
	['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
	['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
	['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], 
											  ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],
											  ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],
	['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], 
											  ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],
											  ['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'],['LDA','PBE'], ['LDA','PBE'],['LDA','PBE'],['LDA','PBE']
	]							   

	pseudo_type_SIESTA = [  ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], 
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],																																																								 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],																																			 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],																																			 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], 
				     					 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
										 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], 
										 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
										 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R']
	]							   

	return pt_atom_num, Atom_Name, pseudo_type_SIESTA, Atom_BG_color

def Pseudo_DB():

	pt_atom_num, Atom_Name, pseudo_type, Atom_BG_color = atom_info_html_table()

	html_string= """ <title>Pseudo potential database</title>
<head> 

<script src="/static/js/three/libs/Three_R86.js" type="text/javascript"></script>
<script src="/static/js/three/libs/TrackballControls_R68.js" type="text/javascript"></script>

<script src="static/js/MarchingCubes/MarchingCubes.js"></script>

<script>
var atomAnalyzerNamespace;
var disabled=false;

var Old_target_x, Old_target_y, Old_target_z;
var Old_camera_x, Old_camera_y, Old_camera_z;

var target_x, target_y, target_z;
var camera_x, camera_y, camera_z;

var imagePath = parent.atomTransitorAnalyzerImagePath;

var scene_width = 370;
var scene_height= 370;

var scene = new THREE.Scene();
var renderer = new THREE.WebGLRenderer();


var raycaster = new THREE.Raycaster();

var Atoms = new THREE.Object3D();
var Link_Group = new THREE.Object3D();
var aBoxLine = new THREE.Object3D();

var F_ini_Box=0;

var BoxLines ;

var BoxLineGroup = new THREE.Object3D();

var BoxGroup = new THREE.Object3D();

var camera;
// = new THREE.PerspectiveCamera(45, scene_width / scene_height, 0.01, 10000);
//var camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 1, 10000 );

var ambientLight = new THREE.AmbientLight(0x1c1c1c); scene.add(ambientLight);

var directionalLight = new THREE.DirectionalLight(0xffffff,1);
//directionalLight.position.set(1, 0, 1).normalize();
scene.add(directionalLight);

var directionalLight2 = new THREE.DirectionalLight(0xffffff,1);
//directionalLight.position.set(1, 0, 1).normalize();
scene.add(directionalLight2);

var mouse = new THREE.Vector2(), INTERSECTED;

var Nx, Ny, Nz;

var xmx, ymx, zmx;
var xmn, ymn, zmn;
var fdf_text;

var N01_P;
var N01_R;


var Dxx, Dyy, Dzz;

var DD_xyz;

var dx, dy ,dz;
var dx_g, dy_g, dz_g;

var Real_N_Atoms ;
var N_Atoms ;
var N_Species, T_N_Species ;

var Mx_Species_Number;

var Lattice_C ;
var unit ;
var R_Bohr = 0.53;

var um = 1e-6;

var N_t=1 ;

var P_Atoms ;
var ID_Atom ;

var N_links;
var Links= new Array(1) ;

var atom_Geo = new Array();
var atom_Mat = new Array();
var atom = new Array();
var Link = new Array();

var Addition_atom;
var ON_save_P=0;


var t_i=0;
var Vertices = {};

var Max_N_Species =50;

var Species_DelOrder = new Array(Max_N_Species);
var Species_Order = new Array(Max_N_Species);
var Species_Number= new Array(Max_N_Species);
var Species_Name = new Array(Max_N_Species);

var iid;
var insec_name;

var exe_mod;

var Extend = new Array(3);
Extend[0]=1;
Extend[1]=1;
Extend[2]=1;

var LatticeVector = new Array(3);

var nonz_i = new Array();

var text_x ;
var text_y ;
var text_z ;


var text_replace_Atom;

var struc_string;

var camera_Vscale=0.7;

var F_exi_Dirtext=1;

var	pseudo_type = [  ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], 
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],																																																								 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],																																			 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],																																			 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], 
					              						  ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
								                 		  ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
	                 ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], 
					               						  ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],
								                 		  ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'], ['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R'],['LDA_N', 'LDA_R', 'GGA_N', 'GGA_R']
	                 ]							   



var Atom_Name=["Removed",
			   "H" , "He",
			   "Li","Be", "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
			   "Na","Mg", "Al","Si","P" ,"S" ,"Cl","Ar",
			   "K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
			   "Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
			   "Cs","Ba",
						 "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
							  "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
			   "Fr","Ra",
						 "Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
							  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"  ];


var Atom_color = [ '000000',
				   'ff4444', '777777',
				   'ff4422','dd2200', 'ffff44','44ff44','ff8800','4444dd','4444ff','666666',
				   'ff2244','dd3300', 'ffff22','22ff22','ffcc00','4400dd','2222ff','555555',
				   'ff2222','dd3300','aa6600','ff9900','ffcc00','ff0000','ff0000','ff0000','ffcc00','ff0000','ff0000','ff0000','ffff00','00ff00','ffcc00','0044dd','0000ff','444444',
				   'ff0000','dd3300','aa6600','ff9900','ffcc00','ff0000','ff0000','ff0000','ffcc00','ff0000','ff0000','ff0000','dddd00','00dd00','ffcc00','2222dd','0000dd','333333',
				   'dd0000','dd0000',
									 'ff0000','dd3300','aa6600','ff9900','ffcc00','ff0000','ff0000','ff0000','ffcc00','ff0000','ff0000','ff0000', 'ff6600','ff9900','ffcc00',
											  'ff0000','ff3300','ff6600','ff9900','ffcc00','ff0000','ff0000','ff0000','ffcc00','bbbb00','00bb00','ff0000', '2200dd','0000bb','222222',
				   'bb0000','dd0000',
									 'ff0000','dd3300','aa6600','ff9900','ffcc00','ff0000','ff0000','ff0000','ffcc00','ff0000','ff0000','ff0000', 'ff6600','ff9900','ffcc00',
											  'ff0000','ff3300','ff6600','ff9900','ffcc00','ff0000','ff0000','ff0000','ffcc00','999900','009900','ff0000', '0022dd','000099','111111'
				   ];

var Atom_R = [ ];

for (var i=1; i<=118; i++) Atom_R[i]=i*0.01+0.15;

var Is_N_spe = new Array(119);
for ( var i=0; i<119; i++ ) Is_N_spe[i]=0;

var Addition_Atom_i, Gra_Atom_i;

var xi, yi, zi;

var x_cam=10.5, y_cam=-10.5, z_cam=10.5;

var Grad = [];
var Rhos = {};

var ang = Math.PI*15/180;

var AtomicCoordinatesFormat ;

var canvas ;

var trackballControl ;

var container_elem;

var plateViewers = new Array(3);
var arrow_plates = new Array(3);

var effect, resolution;

var EP_i, Vr_i;

var V1D_r, V1D_V ; 

var orb_i;

var N_orb, N_lq;
var Text_Info;

var fmin, fmax;

var orb_list2 ;

var F_orb = 0;

var Plate = function() {
	this.mesh = function( mesh ) {
		if ( mesh ) this.Mesh = mesh;
		else return this.Mesh;
	};
	this.arrow = function( arrow ) {
		if ( arrow ) this.Arrow = arrow;
		else return this.Arrow;
	};
	this.geometry = function( geometry ) {
		if ( geometry ) this.Geometry = geometry;
		else return this.Geometry;
	};
	this.material = function( material ) {
		if ( material ) this.Material = material;
		else return this.Material;
	};
	this.vmode = function( mode ) {
		if ( arguments.length === 1 ) this.VMode = mode;
		else return this.VMode;
	};
	this.visibleArrow = function( visible ) {
		var arrow = this.arrow();
		if ( arguments.length === 1 )
		{
			arrow.visible = visible;
			if ( visible ) {
				arrow.geometry.verticesNeedUpdate = true;
			}
		}
		else return arrow.visible;
	};
};

var Canvas = function( container )
{
	var ti=0;

	// From --- middle -----

	this.plates = function( plates ) {
		if ( plates ) this.Plates = plates;
		else return this.Plates;
	};
	this.addPlate = function( id, plate ) {
		var plates = this.plates();
		if ( !plates ) {
			plates = {};
			this.plates(plates);
		}

		plates[id] = plate;
	};
	this.getPlate = function( id ) {
		var plates = this.plates();
		if ( !plates ) return false;

		return plates[id];
	};
	this.setPlateArrow = function( id, arrow ) {
		var plate = this.getPlate(id);
		plate.arrow( plate );
	};
	this.setPlateMesh = function( id, mesh ) {
		var plate = this.getPlate(id);
		plate.mesh( mesh );
	};
	this.setPlateGeometry = function( id, geometry ) {
		var plate = this.getPlate(id);

		plate.geometry( geometry );
	};
	this.setPlateMaterial = function( id, material ) {
		var plate = this.getPlate(id);

		plate.material( material );
	};
	this.getPlateArrow = function( id ) {
		var plate = this.getPlate(id);

		return plate.arrow();
	};
	this.getPlateMesh= function( id ) {
		var plate = this.getPlate(id);

		return plate.mesh();

	};
	this.getPlateGeometry = function( id ) {
		var plate = this.getPlate(id);

		return plate.geometry();
	};
	this.getPlateMaterial = function( id ) {
		var plate = this.getPlate(id);

		return plate.material();
	};

	container.appendChild( renderer.domElement );
	
};

var animate = function() {
	
	Old_target_x = trackballControl.target.x;	Old_target_y = trackballControl.target.y;	Old_target_z = trackballControl.target.z;
	Old_camera_x = camera.position.x;			Old_camera_y = camera.position.y;			Old_camera_z = camera.position.z;
		
	trackballControl.update();

	if(Old_camera_x!=camera.position.x || Old_camera_y!=camera.position.y|| Old_camera_z!=camera.position.z || trackballControl.target.x !=Old_target_x || trackballControl.target.y !=Old_target_y || trackballControl.target.z !=Old_target_z )
	{		
		camera_x = camera.position.x		; camera_y = camera.position.y		; camera_z = camera.position.z;
		target_x = trackballControl.target.x; target_y = trackballControl.target.y; target_z = trackballControl.target.z;
		
	
		if(exe_mod==0) document.getElementById("Views").value = getViews();
	 
	   camera.up = new THREE.Vector3(0.0, 0.0, 1.0);	  
	}
	
	renderer.render(scene, camera );
	requestAnimationFrame(animate);

	};

var string_to_DFT="";
function select_pseudo(a_name, a_color)
{

	var p_name = "pseudo_"+ a_name;				  
	var p_elem = document.getElementById(p_name);	
	
	var wrap_pseu = document.getElementById("pseu_sel");	  

	if (p_elem.selectedIndex !=0)
	{
		create_pseu_tex(a_name, p_name, p_elem, a_color);		
	}
	else // if sel=0, initialize 
	{
		var del_pseu = document.getElementById("sel" + p_name);	 
		del_pseu.remove();

		p_elem.style.fontSize ="13px";
		p_elem.style.color ="#000000";
		p_elem.style.backgroundColor ="#ffffff";
		p_elem.removeAttribute("name");
		document.getElementById('id_curr_PP_info').value = "Pseudo potential :" ;	
		View_TextInfo(1);
		View_Vr(1);
		removeAllObjects();
		View_CS(1);
	} 

	if (wrap_pseu.hasChildNodes())
	{ 
		string_to_DFT="";
		var children = wrap_pseu.childNodes; 
		//console.log('children.length=', children.length);
	  for(var i=2; i<children.length; i++)
	  { 
		  string_to_DFT += ":" + children[i].value ;	
		 // console.log('children[i].value=', children[i].value);	  
	  } 
		  string_to_DFT += "@";
	  for(var i=2; i<children.length; i++)
	  { 
		  string_to_DFT += ":" + children[i].style.backgroundColor ;		  
	  } 

	}


	show_info(p_elem, 0);

}

function create_pseu_tex(a_name, p_name, p_elem, a_color)
{
		p_elem.style.fontSize ="13px";
		p_elem.style.color ="#ffffff";
		p_elem.style.backgroundColor ="#555555";

		
		var wrap_pseu = document.getElementById("pseu_sel");											 

		var new_pseu = document.createElement("input");
		new_pseu.type="text";
		new_pseu.value = a_name + ' ' + p_elem.options[p_elem.selectedIndex].text  ;
		wrap_pseu.appendChild(new_pseu);		
		
		new_pseu.name = "pseudo";

		new_pseu.id ="sel" + p_name;

		new_pseu.draggable = "true";		
		new_pseu.class = "dropzone" ;	   
        new_pseu.ondragstart = "event.dataTransfer.setData('text/plain',null)";

        new_pseu.style.cursor = "pointer";
		new_pseu.style.textAlign = "center";
		new_pseu.style.color = "#ffffff";
		new_pseu.style.backgroundColor = "#"+a_color;
		new_pseu.style.border = "3px solid #ccccff";
						
		new_pseu.size = "7"
		new_pseu.setAttribute("readonly","");

		//write to curr textbox
		document.getElementById('id_curr_PP_info').value = "Selected pseudo potential : " + new_pseu.value ;	



}


var dragged;

/* 드래그 가능한 대상에서 이벤트 발생 */
document.addEventListener("drag", function(event) {

}, false);

// 요소를 반투명하게 함
document.addEventListener("dragstart", function(event) { dragged = event.target; event.target.style.opacity = .5; }, false);

// 투명도를 리셋
document.addEventListener("dragend", function(event) { event.target.style.opacity = ""; }, false);

/* 드롭 대상에서 이벤트 발생 */
document.addEventListener("dragover", function(event) { event.preventDefault(); }, false);

// 요소를 드롭하려는 대상 위로 드래그했을 때 대상의 배경색 변경
document.addEventListener("dragenter", function(event) { if (event.target.className == "dropzone") { event.target.style.background = "purple"; } }, false);

// 요소를 드래그하여 드롭하려던 대상으로부터 벗어났을 때 배경색 리셋
document.addEventListener("dragleave", function(event) { if (event.target.className == "dropzone") { event.target.style.background = "";  } }, false);

// 기본 액션을 막음 (링크 열기같은 것들)
// 드래그한 요소를 드롭 대상으로 이동
document.addEventListener("drop", function(event) 
{ 
    event.preventDefault(); 
    target_p = event.target
  
    if (event.target.class == "dropzone") 
    {       

        var dummy = target_p.value;  target_p.value = dragged.value; dragged.value = dummy;
        var dummy = target_p.style.backgroundColor;  
		target_p.style.backgroundColor = dragged.style.backgroundColor; 
		dragged.style.backgroundColor  = dummy;

        var dummy = target_p.id; target_p.id = dragged.id; dragged.id = dummy;        

      

	  	var wrap_pseu = document.getElementById("pseu_sel");	  


		if (wrap_pseu.hasChildNodes())
		{ 
			string_to_DFT="";
			var children = wrap_pseu.childNodes; 
			console.log('children.length=', children.length);
		for(var i=2; i<children.length; i++)
		{ 
			string_to_DFT += ":" + children[i].value ;	
			console.log('children[i].value=', children[i].value);	  
		} 
			string_to_DFT += "@";
		for(var i=2; i<children.length; i++)
		{ 
			string_to_DFT += ":" + children[i].style.backgroundColor ;		  
		} 

		}

       
    }
  
}, false);

function send_to_DFT() {
	localStorage.setItem('Load_PP_info', string_to_DFT);
	//event.dataTransfer.setData("text/plain", string_to_DFT);
}

function loadEPData( data ) {
	
	exe_mod=1;

	container_elem = document.getElementById("container");
	canvas = new Canvas( container_elem );

	camera = new THREE.PerspectiveCamera(45, scene_width / scene_height, 0.01, 10000);

	camera.up = new THREE.Vector3(0, 0, 1);
	
	
	renderer.clear(true,true)

	renderer.setClearColor( 0xffffff, 1 );
	renderer.setSize(scene_width, scene_height);
	renderer.shadowMapEnabled = true;
	renderer.sortObjects = true;
	
	removeAllObjects();

	var fin_i=0;	

	var lines = data.split('\\n');
	
	ON_save_P=1;	

	xmx= Number(-1e+100);
	ymx= Number(-1e+100);
	zmx= Number(-1e+100);
	xmn= Number(1e+100);
	ymn= Number(1e+100);
	zmn= Number(1e+100);

	
	for ( var index=0; index<lines.length; index++ ) {
		
		var line = lines[index].trim();		

		if ( !line ) continue;		

		if ( line[0] === '@') {
			var line = line.replace('@', '');
			var keyVal = line.split(' ');
			var start = index+1;
			var end = start + Number(keyVal[1]);
			index = end-1;
			switch( keyVal[0].trim()) {				
				case 'EP_data'  : parseEPData( lines.slice(start, end)); break;									
				case 'V_r'	    : if(F_orb==0) parse_Vr( lines.slice(start, end)); break;	
				case 'orbs'	    : if(F_orb==0) parse_orbs( lines.slice(start, end)); break;
				case 'textinfo' : if(F_orb==0) parseTextInfo( lines.slice(start, end)); break;				
				case 'params'   : if(F_orb==0) parseParams( lines.slice(start, end)); break;
							
				default: alert( 'File type mismatch......');
			}
		}
		
	}

	

	View_Vr(0);
	
	View_TextInfo(0);

	xi=1; yi=1; zi=1;	

	text_x = document.getElementById("position_x");
	text_y = document.getElementById("position_y");
	text_z = document.getElementById("position_z");
	
	directionalLight.position.set(3*xmx, 3*ymx, 3*zmx).normalize();
	directionalLight2.position.set(-3*xmx, -3*ymx, -3*zmx).normalize();

	trackballControl = new THREE.TrackballControls(camera, renderer.domElement);

	trackballControl.rotateSpeed = DD_xyz*0.05;
	console.log('trackballControl.rotateSpeed=', trackballControl.rotateSpeed);
	trackballControl.zoomSpeed = 1.0;
	trackballControl.panSpeed = 1.0;

	camera.rotation.z=Math.PI*0.5;

	camera.position.set(camera_Vscale*(Dyy+Dzz)+Dxx*0.5, camera_Vscale*(Dzz+Dxx)+Dyy*0.5, camera_Vscale*(Dxx+Dyy)+Dzz*0.5);
	trackballControl.target = new THREE.Vector3((Dxx)*0.5, (Dyy)*0.5, (Dzz)*0.5);
	
	var plate = new Plate();
	arrow_plates[0] = new THREE.LineSegments( new THREE.Geometry(), new THREE.LineBasicMaterial( { color: 0x000000, opacity: 1, linewidth: 1 } ) );
	arrow_plates[0].position.z=0;
	arrow_plates[0].visible = false;

	plate.arrow( arrow_plates[0] );
	plate.vmode(0);
	canvas.addPlate( 'xy', plate );
	scene.add( arrow_plates[0] );

	plate = new Plate();
	arrow_plates[1] = new THREE.LineSegments( new THREE.Geometry(), new THREE.LineBasicMaterial( { color: 0x000000, opacity: 1, linewidth: 1 } ) );
	arrow_plates[1].position.x=0;
	arrow_plates[1].visible = false;
	plate.arrow( arrow_plates[1] );
	plate.vmode(0);
	canvas.addPlate( 'yz', plate );
	scene.add( arrow_plates[1] );

	plate = new Plate();
	arrow_plates[2] = new THREE.LineSegments( new THREE.Geometry(), new THREE.LineBasicMaterial( { color: 0x000000, opacity: 1, linewidth: 1 } ) );
	arrow_plates[2].position.y=0;
	arrow_plates[2].visible = false;
	plate.arrow( arrow_plates[2] );
	plate.vmode(0);
	canvas.addPlate( 'zx', plate );
	scene.add( arrow_plates[2] );



	var plateGeometry = new THREE.Geometry();
	var plateMaterial = new THREE.MeshBasicMaterial({ vertexColors: THREE.VertexColors, transparent:true , opacity: 0.5, side:THREE.DoubleSide });
	plateViewers[0] = new THREE.Mesh( plateGeometry, plateMaterial );
	plateViewers[0].position.z = 0;
	plateViewers[0].visible=false;
	canvas.setPlateGeometry('xy', plateGeometry);
	canvas.setPlateMaterial('xy', plateMaterial);
	canvas.setPlateMesh('xy', plateViewers[0]);
	scene.add( plateViewers[0] );

	plateGeometry = new THREE.Geometry();
	plateMaterial = new THREE.MeshBasicMaterial({ vertexColors: THREE.VertexColors, transparent:true , opacity: 0.5, side:THREE.DoubleSide });
	plateViewers[1] = new THREE.Mesh( plateGeometry, plateMaterial );
	plateViewers[1].position.x = 0;
	plateViewers[1].visible=false;
	canvas.setPlateGeometry('yz', plateGeometry);
	canvas.setPlateMaterial('yz', plateMaterial);
	canvas.setPlateMesh('yz', plateViewers[1]);
	scene.add( plateViewers[1] );

	plateGeometry = new THREE.Geometry();
	plateMaterial = new THREE.MeshBasicMaterial({ vertexColors: THREE.VertexColors, transparent:true , opacity: 0.5, side:THREE.DoubleSide });
	plateViewers[2] = new THREE.Mesh( plateGeometry, plateMaterial );
	plateViewers[2].position.y = 0;
	plateViewers[2].visible=false;
	canvas.setPlateGeometry('zx', plateGeometry);
	canvas.setPlateMaterial('zx', plateMaterial);
	canvas.setPlateMesh('zx', plateViewers[2]);
	scene.add( plateViewers[2] );
	
	if (ON_save_P==1)
	{

	setViewerShape('xy');
	setViewerShape('yz');
	setViewerShape('zx');

	setViewerFaces();

	setArrows('xy');
	setArrows('yz');
	setArrows('zx');
		
	}

	
	
	drawBox(); //<----don't delete 

	//document.removeEventListener( 'click', onDocumentMouseMove, false );
	//document.addEventListener( 'mousemove', onDocumentMouseMove, false );
	//document.addEventListener( 'mouseover', onDocumentMouseOver, false );
	//window.addEventListener( 'resize', onWindowResize, false );
	
	

	//	dropZone1.addEventListener('dragstart', handleDragStart, false);




	
	MC_effect();
		
			
	animate();
	
}

function MC_effect()
{
	var color_MC= new THREE.Color( 0.9, 0.9, 0.9 );
	var materials = new THREE.MeshStandardMaterial({opacity: 0.4, color: color_MC, transparent:true, side:THREE.DoubleSide}) ;
	var resolution = 40;

	var Equi_v = document.getElementById('MC_Range').value;

	document.getElementById('MC_Range_text').value = Equi_v;

	if (typeof(effect)!="undefined") {effect.reset();}

	scene.remove( effect );

	effect = new THREE.MarchingCubes( resolution, materials, false, true );
	effect.position.set( xmx/2, ymx/2, zmx/2 );
	effect.scale.set(  xmx*0.5,  ymx*0.5,  zmx*0.5 );

	scene.add( effect );

	effect.reset();

	var subtract, strength;

		subtract = 5.0;			
		strength = 0.05;

	var vertices = Vertices[0];	
	
	for(var i=0; i<Nx; i++)	{
	for(var j=0; j<Ny; j++)	{
	for(var k=0; k<Nz; k++)	{
				
		var ord_xyz = k + j*Nz + i*Nz*Ny;
		if (vertices[ord_xyz].p<Equi_v)
		{			
			effect.addBall(vertices[ord_xyz].x/xmx, vertices[ord_xyz].y/ymx, vertices[ord_xyz].z/zmx, strength, subtract);			
		}
	}}}

}

function removeAllObjects()
{
	
	aBoxLine = new THREE.Object3D();

	scene.remove( BoxLineGroup ); BoxLineGroup = new THREE.Object3D();

	scene.remove(BoxGroup); BoxGroup = new THREE.Object3D();
	for (var i=0; i<3;i++){ scene.remove(plateViewers[i]); scene.remove(arrow_plates[i]);}
	if (typeof(effect)!="undefined") {effect.reset();}
	scene.remove(effect);

}
function ON_MouseMove(){ document.addEventListener( 'mousemove', onDocumentMouseMove, false );}

function OFF_MouseMove(){ document.removeEventListener( 'mousemove', onDocumentMouseMove, false ); }

function parseTextInfo ( dataLines ) {

Text_Info  = dataLines[0] +"\\n";
	Text_Info += dataLines[1] +"\\n";

	var psf_orb = dataLines[2].split('/');
	
	for (var i=0; i< psf_orb.length-1; i++)
	{
		Text_Info += psf_orb[i] +"\\n" ;
		
	}

	Text_Info += dataLines[3].trim() +"\\n";


}


function parse_Vr ( dataLines ) {

	
	V1D_r[Vr_i] = new Array(0); 
	V1D_V[Vr_i] = new Array(0); 
		
	for ( var i=1; i<dataLines.length ; i++ ) 
	{
		var line   = dataLines[i];
		var values = line.split(' ');
		
		V1D_r[Vr_i].push(values[0].trim());
		V1D_V[Vr_i].push(values[1].trim());
	
	};		
		

	Vr_i++;

	
};

var orb_name;

function parse_orbs ( dataLines ) {
	

	var orb_list = document.getElementById('id_orb_list');
	
	while ( orb_list.hasChildNodes() ) { orb_list.removeChild( orb_list.firstChild ); }

	var orb_opt = new Array(N_orb);

	orb_name = new Array(N_orb);
			
	for ( var i=0; i<dataLines.length ; i++ ) 
	{
		var line   = dataLines[i];
		var values = line.split('orb: ');
		
		orb_name[i] = values[1];

		orb_opt[i] = document.createElement("option");
		orb_opt[i].value = String(i+1) ;
		orb_opt[i].text = orb_name[i];
	
		orb_list.add(orb_opt[i], null);
	
	};		

};


function parseParams ( dataLines ) {	
	

		for ( var i=0; i<dataLines.length; i++ ) {
			var line = dataLines[i];

			var values = line.split('=');

			switch( values[0].trim() ) {
				case 'Nx' :Nx =Number(values[1]); break;
				case 'Ny' :Ny =Number(values[1]); break;
				case 'Nz' :Nz =Number(values[1]); break;
				case 'N_nod' :N_nod =Number(values[1]); break;
				case 'N_t' :N_t =Number(values[1]); break;
				case 'dx' :dx_g =Number(values[1]); break;
				case 'dy' :dy_g =Number(values[1]); break;
				case 'dz' :dz_g =Number(values[1]); break;
				case 'N_orb' :N_orb =Number(values[1]); break;		
				case 'N_lq' :N_lq =Number(values[1]); break;		
				default: alert('Un-recognizable parameter: '+values[0].trim()); return;
			}

			
		};		

		

		
		V1D_r = new Array(N_lq);
		V1D_V = new Array(N_lq);


		Vr_i = 0;

		document.getElementById('PxyRange').max = Nz ;
		document.getElementById('PyzRange').max = Nx ;
		document.getElementById('PzxRange').max = Ny ;


	};

			

	function parseEPData ( dataLines ) {

//--------------- to get V field  --------------

		EP_i++;

		orb_dummy = dataLines[0].split(' ');
		orb_dummy = orb_dummy[0].split(':');

		if (orb_i != orb_dummy[1]) return;

		if (ON_save_P===0) return;

		var time = 0;
		var xi_L, yi_L, zi_L, md;
		var c1=0.25, c2=0.50, c3=0.75;
		var R_col, G_col, B_col;
		var func= new Array(N_nod);

		Vertices[time] = [];
		Grad = [];	

		var Nzy = Nz*Ny;

		var fminmax = dataLines[0].split('#');

		fmin = fminmax[1];
		fmax = fminmax[2];
	

		for ( var i=0; i<dataLines.length-1; i++ ) {
			var line = dataLines[i+1];
			

			xi_L = parseInt(i/Nzy); md = parseInt(i%Nzy) ;
			yi_L = parseInt(md/Nz);
			zi_L = parseInt(md%Nz);

			func[i]=parseFloat(line);

			if (func[i]>c3 && func[i]<=1.0) { R_col= 1.0 ; G_col= (1.0-func[i])/(1.0-c3) ; B_col= 0.0 ; }
	   else if (func[i]>c2 && func[i]<=c3 ) { R_col= (func[i]-c2)/(c3-c2) ; G_col= 1.0 ; B_col= 0.0 ; }
	   else if (func[i]>c1 && func[i]<=c2 ) { R_col= 0.0 ; G_col= 1.0 ; B_col= (c2-func[i])/(c2-c1) ; }
	   else if (func[i]>=0 && func[i]<=c1 ) { R_col= 0.0 ; G_col= (func[i]-0.0)/(c1-0.0) ; B_col= 1.0 ; }

			var vertex = {};
			vertex.x = xi_L*dx_g;
			vertex.y = yi_L*dy_g;
			vertex.z = zi_L*dz_g;
			vertex.cr = R_col;
			vertex.cg = G_col;
			vertex.cb = B_col;
			vertex.p = func[i];

			if (vertex.x>xmx) xmx = vertex.x;
			if (vertex.y>ymx) ymx = vertex.y;
			if (vertex.z>zmx) zmx = vertex.z;

			if (vertex.x<xmn) xmn = vertex.x;			
			if (vertex.y<ymn) ymn = vertex.y;			
			if (vertex.z<zmn) zmn = vertex.z;

			Vertices[time].push( vertex );

		}

		
		Dxx = xmx-xmn; Dyy = ymx-ymn; Dzz = zmx-zmn;

		DD_xyz = Math.sqrt(Dxx*Dxx + Dyy*Dyy + Dzz*Dzz);

		//--------------- to get E field  --------------

		var ord = new Array(Nx);

		for (var i=0; i<Nx; i++) ord[i] = new Array(Ny);
		for (var i=0; i<Nx; i++)
		for (var j=0; j<Ny; j++) ord[i][j] = new Array(Nz);

		var grad_fx = new Array(N_nod);
		var grad_fy = new Array(N_nod);
		var grad_fz = new Array(N_nod);

		var grad_fx_mx=-1e+100, grad_fy_mx=-1e+100, grad_fz_mx=-1e+100;
		var grad_fx_mn= 1e+100, grad_fy_mn= 1e+100, grad_fz_mn= 1e+100;

		for (var i=0;i<N_nod;i++)
		{
			var dummy;
			xi_L = parseInt(i / Nzy) ; dummy = i % Nzy;
			yi_L = parseInt(dummy / Nz) ;
			zi_L = dummy % Nz ;
			ord[xi_L][yi_L][zi_L] = zi_L + Nz*yi_L + Nz*Ny*xi_L ;
		}

		var Norm_dx = dx_g/(xmx-xmn) ;
		var Norm_dy = dy_g/(ymx-ymn) ;
		var Norm_dz = dz_g/(zmx-zmn) ;

		for (var i=0;i<N_nod;i++)
		{
			var dummy;

			xi_L = parseInt(i / Nzy) ; dummy = i % Nzy;
			yi_L = parseInt(dummy / Nz) ;
			zi_L = dummy % Nz ;

		 	if (xi_L==0 || xi_L==Nx-1) grad_fx[i] = 0; else grad_fx[i] = ( func[ord[xi_L+1][yi_L][zi_L]] - func[ord[xi_L-1][yi_L][zi_L]] ) / dx_g;
			if (yi_L==0 || yi_L==Ny-1) grad_fy[i] = 0; else grad_fy[i] = ( func[ord[xi_L][yi_L+1][zi_L]] - func[ord[xi_L][yi_L-1][zi_L]] ) / dy_g;
			if (zi_L==0 || zi_L==Nz-1) grad_fz[i] = 0; else grad_fz[i] = ( func[ord[xi_L][yi_L][zi_L+1]] - func[ord[xi_L][yi_L][zi_L-1]] ) / dz_g;

			if (grad_fx[i]> grad_fx_mx) grad_fx_mx = grad_fx[i]; if (grad_fx[i]< grad_fx_mn) grad_fx_mn = grad_fx[i];
			if (grad_fy[i]> grad_fy_mx) grad_fy_mx = grad_fy[i]; if (grad_fy[i]< grad_fy_mn) grad_fy_mn = grad_fy[i];
			if (grad_fz[i]> grad_fz_mx) grad_fz_mx = grad_fz[i]; if (grad_fz[i]< grad_fz_mn) grad_fz_mn = grad_fz[i];

		}

		var diffx_mx = grad_fx_mx - grad_fx_mn;
		var diffy_mx = grad_fy_mx - grad_fy_mn;
		var diffz_mx = grad_fz_mx - grad_fz_mn;

		var grad_lim=16.0;

		for (i=0;i<N_nod;i++)
		{

			var dummy;

			var Norm_grad_x, Norm_grad_y, Norm_grad_z;

			xi_L = parseInt(i / Nzy) ; dummy = i % Nzy;
			yi_L = parseInt(dummy / Nz) ;
			zi_L = dummy % Nz ;

			Norm_grad_x = -grad_fx[i]/diffx_mx ; if (Math.abs(Norm_grad_x)>grad_lim*Norm_dx) Norm_grad_x = Norm_grad_x/Math.abs(Norm_grad_x)*grad_lim*Norm_dx ;
			Norm_grad_y = -grad_fy[i]/diffy_mx ; if (Math.abs(Norm_grad_y)>grad_lim*Norm_dy) Norm_grad_y = Norm_grad_y/Math.abs(Norm_grad_y)*grad_lim*Norm_dy ;
			Norm_grad_z = -grad_fz[i]/diffz_mx ; if (Math.abs(Norm_grad_z)>grad_lim*Norm_dz) Norm_grad_z = Norm_grad_z/Math.abs(Norm_grad_z)*grad_lim*Norm_dz ;


			var xyz = {
				x: Norm_grad_x,
				y: Norm_grad_y,
				z: Norm_grad_z
			};
			Grad.push( xyz );
		}

	};

		

//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

//----------------------------------------
	var V_2D ;
		
	function setViewerShape ( plateId )
	{

		
		var outMax, inMax, order, xyzp;
		var Nzy = Nz*Ny;
		switch( plateId ) {
			case 'xy':
				outMax = Nx;
				inMax = Ny;
				order = function(i, j){ return (zi-1) + Nz*(j-1) + Nzy*(i-1); };
				xyzp = function( i, j ){ return (j-1 ) + Ny*(i-1); };
				break;
			case 'yz':
				outMax = Ny;
				inMax = Nz;
				order = function(j, k){ return (k-1) + Nz*(j-1) + Nzy*(xi-1); };
				xyzp = function( j, k ){ return (k-1) + Nz*(j-1); };
				break;
			case 'zx':
				outMax = Nz;
				inMax = Nx;
				order = function(k, i){ return (k-1) + Nz*(yi-1) + Nzy*(i-1); };
				xyzp = function( k, i ){ return (i-1) + Nx*(k-1); };
				break;
		}

		var plate = canvas.getPlate( plateId );
		var geometry = plate.geometry();		

		var vertices = Vertices[0];


		V_2D = new Array(outMax);
		
		for (var i = 0; i < outMax; i++) { V_2D[i] = new Array(inMax); }
		for (var i=1; i<=outMax ; i++) {
			for (var j=1; j<=inMax ; j++) {
				var ord_data = order( i, j );
				var ord_XYP  = xyzp( i, j );

				geometry.vertices[ord_XYP]=	new THREE.Vector3(
									  vertices[ord_data].x + (plateId=="yz"? 1:0)* vertices[ord_data].p*plate.vmode(),
									  vertices[ord_data].y - (plateId=="zx"? 1:0)* vertices[ord_data].p*plate.vmode(),
									  vertices[ord_data].z + (plateId=="xy"? 1:0)* vertices[ord_data].p*plate.vmode() ) ;
			  
				 V_2D[i-1][j-1] = vertices[ord_data].p;
			}
		}




	};


//---------------------------------------------------------------------------------------------------------
	function setViewerFaces () {
		var geoXY = canvas.getPlateGeometry('xy');
		var geoYZ = canvas.getPlateGeometry('yz');
		var geoZX = canvas.getPlateGeometry('zx');

		for (var i=1; i<=Nx-1; i++) {
			for (var j=1; j<=Ny-1; j++) {
				var ord_XYP = (j-1) + Ny*(i-1);

				geoXY.faces.push( new THREE.Face3(ord_XYP , ord_XYP+1, ord_XYP + Ny ) );
				geoXY.faces.push( new THREE.Face3(ord_XYP + Ny, ord_XYP+1, ord_XYP + Ny+1) );

			}
		}

			//--------------------------------------- YZ �� ---------------------------

		for (var j=1; j<=Ny-1; j++) {
			for (var k=1; k<=Nz-1; k++) {
				var ord_YZP = (k-1) + Nz*(j-1);

				geoYZ.faces.push( new THREE.Face3(ord_YZP , ord_YZP + 1 , ord_YZP+Nz ) );
				geoYZ.faces.push( new THREE.Face3(ord_YZP + Nz, ord_YZP + 1, ord_YZP+Nz+1) );

			}
		}

			//--------------------------------------- ZX �� ---------------------------

		for (var k=1; k<=Nz-1; k++) {
			for (var i=1; i<=Nx-1; i++) {
				var ord_ZXP = (i-1) + Nx*(k-1);

				geoZX.faces.push( new THREE.Face3(ord_ZXP , ord_ZXP+1, ord_ZXP + Nx ) );
				geoZX.faces.push( new THREE.Face3(ord_ZXP + Nx, ord_ZXP+1, ord_ZXP + Nx+1) );
			}
		}
	};
function setArrows( plateId )
{
		var plate = canvas.getPlate( plateId );
		var i_arr=0;
		var arrow_L = 0.5;
		var s_L=0.4;
		var expan=2;

		var arrow = canvas.getPlateArrow( plateId );
		var geometry = arrow.geometry;

		var outMax, inMax, order;

		var Nzy = Nz*Ny;

		var vertices = Vertices[0];

		switch( plateId ) {
			case 'xy':
				outMax = Nx;
				inMax = Ny;
				order = function(i, j){ return (zi-1) + Nz*j + Nzy*i ; };

				for (var i=0; i<outMax; i=i+expan) {
					for (var j=0; j<inMax; j=j+expan) {
						var ord = order(i, j);

						var x0=vertices[ord].x - Grad[ord].x*arrow_L*expan;
						var y0=vertices[ord].y - Grad[ord].y*arrow_L*expan;
						var x1=vertices[ord].x + Grad[ord].x*arrow_L*expan;
						var y1=vertices[ord].y + Grad[ord].y*arrow_L*expan;

						geometry.vertices[i_arr++]= new THREE.Vector3(x0, y0, vertices[ord].z) ;
						geometry.vertices[i_arr++]= new THREE.Vector3(x1, y1, vertices[ord].z) ;

						var x2 = ((x0-x1)*Math.cos(ang)-(y0-y1)*Math.sin(ang))*s_L+x1, y2 = ((x0-x1)*Math.sin(ang)+(y0-y1)*Math.cos(ang))*s_L+y1;

						geometry.vertices[i_arr++]= new THREE.Vector3(x1, y1, vertices[ord].z) ;
						geometry.vertices[i_arr++]= new THREE.Vector3(x2, y2, vertices[ord].z) ;

						var x3 = ((x0-x1)*Math.cos(-ang)-(y0-y1)*Math.sin(-ang))*s_L+x1, y3 = ((x0-x1)*Math.sin(-ang)+(y0-y1)*Math.cos(-ang))*s_L+y1;

						geometry.vertices[i_arr++]= new THREE.Vector3(x1, y1, vertices[ord].z) ;
						geometry.vertices[i_arr++]= new THREE.Vector3(x3, y3, vertices[ord].z) ;
					}
				}

				break;
			case 'yz':
				outMax = Ny;
				inMax = Nz;
				order = function(j, k){ return k + Nz*j + Nzy*(xi-1) ; };

				for (var j=0; j<outMax; j=j+expan)
				{
				  for (var k=0; k<inMax; k=k+expan)
				  {
					  var ord = order(j,k);


						 var y0=vertices[ord].y- Grad[ord].y*arrow_L*expan;
						 var z0=vertices[ord].z- Grad[ord].z*arrow_L*expan;
						 var y1=vertices[ord].y+ Grad[ord].y*arrow_L*expan;
						 var z1=vertices[ord].z+ Grad[ord].z*arrow_L*expan;

						 geometry.vertices[i_arr++]= new THREE.Vector3(vertices[ord].x, y0, z0) ;
						 geometry.vertices[i_arr++]= new THREE.Vector3(vertices[ord].x, y1, z1) ;

						 var y2 = ((y0-y1)*Math.cos(ang)-(z0-z1)*Math.sin(ang))*s_L+y1;
						 var z2 = ((y0-y1)*Math.sin(ang)+(z0-z1)*Math.cos(ang))*s_L+z1;

						 geometry.vertices[i_arr++]= new THREE.Vector3(vertices[ord].x, y1, z1) ;
						 geometry.vertices[i_arr++]= new THREE.Vector3(vertices[ord].x, y2, z2) ;

						 var y3 = ((y0-y1)*Math.cos(-ang)-(z0-z1)*Math.sin(-ang))*s_L+y1;
						 var z3 = ((y0-y1)*Math.sin(-ang)+(z0-z1)*Math.cos(-ang))*s_L+z1;

						 geometry.vertices[i_arr++]= new THREE.Vector3(vertices[ord].x, y1, z1) ;
						 geometry.vertices[i_arr++]= new THREE.Vector3(vertices[ord].x, y3, z3) ;

					}
				}
				break;
			case 'zx':
				outMax = Nz;
				inMax = Nx;
				order = function(k, i){ return k + Nz*(yi-1) + Nzy*i ; };

				for (var k=0; k<outMax ; k=k+expan)
				   {
					for (var i=0; i<inMax ; i=i+expan)
					 {
							 var ord = order(k,i);
							 var z0=vertices[ord].z- Grad[ord].z*arrow_L*expan;
							 var x0=vertices[ord].x- Grad[ord].x*arrow_L*expan;
							 var z1=vertices[ord].z+ Grad[ord].z*arrow_L*expan;
							 var x1=vertices[ord].x+ Grad[ord].x*arrow_L*expan;

							 geometry.vertices[i_arr++]= new THREE.Vector3(x0, vertices[ord].y, z0) ;
							 geometry.vertices[i_arr++]= new THREE.Vector3(x1, vertices[ord].y, z1) ;

							 var x2 = ((x0-x1)*Math.cos(ang)-(z0-z1)*Math.sin(ang))*s_L+x1;
							 var z2 = ((x0-x1)*Math.sin(ang)+(z0-z1)*Math.cos(ang))*s_L+z1;

							 geometry.vertices[i_arr++]= new THREE.Vector3(x1, vertices[ord].y, z1) ;
							 geometry.vertices[i_arr++]= new THREE.Vector3(x2, vertices[ord].y, z2) ;

							 var x3 = ((x0-x1)*Math.cos(-ang)-(z0-z1)*Math.sin(-ang))*s_L+x1;
							 var y3 = ((x0-x1)*Math.sin(-ang)+(z0-z1)*Math.cos(-ang))*s_L+z1;

							 geometry.vertices[i_arr++]= new THREE.Vector3(x1, vertices[ord].y, z1) ;
							 geometry.vertices[i_arr++]= new THREE.Vector3(x3, vertices[ord].y, y3) ;
						}
					  }

				break;
		}

		if ( arrow.visible )
			geometry.verticesNeedUpdate = true;
	};

	//------------------------
	function drawBox () {
		
		var BLines = [];		
		
		var Lp_x=xmx, Lp_y=ymx, Lp_z=zmx;
				
		var Lines_p = [[0,0,0],[Lp_x,0,0],[Lp_x,Lp_y,0],[0,Lp_y,0],[0,0,Lp_z],[Lp_x,0,Lp_z],[Lp_x,Lp_y,Lp_z],[0,Lp_y,Lp_z]];

		var LL = [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]];

		var N_Lines = LL.length;

		aBoxLine = new THREE.Object3D();		
		
		for ( var i=0;i<=N_Lines-1; i++)
		{
			BLines[i] = {};
			BLines[i].geo = new THREE.Geometry();

			BLines[i].geo.vertices.push(new THREE.Vector3( Lines_p[LL[i][0]][0], Lines_p[LL[i][0]][1], Lines_p[LL[i][0]][2]));
			BLines[i].geo.vertices.push(new THREE.Vector3( Lines_p[LL[i][1]][0], Lines_p[LL[i][1]][1], Lines_p[LL[i][1]][2]));

			BLines[i].mat = new THREE.LineDashedMaterial( { color: 0x666666, linewidth: 1, scale:1, dashSize: 3, gapSize: 1 } );
			BLines[i].mesh = new THREE.Line( BLines[i].geo, BLines[i].mat );
			

			aBoxLine.add(BLines[i].mesh);
		}
			
										BoxLines	= new Array(Extend[0]);
		for (var i=0; i<Extend[0]; i++) BoxLines[i] = new Array(Extend[1]);
		for (var i=0; i<Extend[0]; i++)
		for (var j=0; j<Extend[1]; j++) BoxLines[i][j] = new Array(Extend[2]);

		scene.remove( BoxLineGroup );
		BoxLineGroup = new THREE.Object3D();

		for (var i=0; i<Extend[0]; i++)
		for (var j=0; j<Extend[1]; j++)
		for (var k=0; k<Extend[2]; k++)
		{
			BoxLines[i][j][k]=aBoxLine.clone();
			BoxLines[i][j][k].position.set(i*Lp_x, j*Lp_y, k*Lp_z);
			BoxLineGroup.add(BoxLines[i][j][k]);
		}

		scene.add( BoxLineGroup );

		// Box is not drawn.  -----------------------------------------

		var aBox = new THREE.Mesh(new THREE.BoxGeometry(Lp_x,Lp_y,Lp_z), new THREE.MeshLambertMaterial({color: 0xffffff, transparent:true, opacity:0.1}));

										var Boxes = new Array(Extend[0]);
			for (var i=0; i<Extend[0]; i++) Boxes[i] = new Array(Extend[1]);
			for (var i=0; i<Extend[0]; i++)
			for (var j=0; j<Extend[1]; j++) Boxes[i][j] = new Array(Extend[2]);

			scene.remove( BoxGroup );
			BoxGroup = new THREE.Object3D();

			for (var i=0; i<Extend[0]; i++)
			for (var j=0; j<Extend[1]; j++)
			for (var k=0; k<Extend[2]; k++)
			{
				Boxes[i][j][k]=aBox.clone();

				Boxes[i][j][k].position.set(i*Lp_x+Lp_x*0.5, j*Lp_y+Lp_y*0.5, k*Lp_z+Lp_z*0.5);

				BoxGroup.add( Boxes[i][j][k] );
			}

	}; // End of drawAxis()
//---------------------------------------------------------------------
function onDocumentMouseMove( event )
{


}
//---------------------------------------------------------------------
function onDocumentMouseOver( event )
{
	event.preventDefault();

}

function onWindowResize() {
	/*
	scene_width = document.body.clientWidth;
	scene_height= document.body.clientHeight;

	console.log('==================================');
	console.log('OnResize: '+document.body.clientWidth + ', ' + document.body.clientHeight);
	console.log('OnResize(): '+scene_width+', '+scene_height);
	
	camera.aspect = scene_width / scene_height;
	camera.updateProjectionMatrix();

	renderer.setSize( scene_width, scene_height );
	*/

}
//-------------------------------

function toggleView( plateId ) {
	setPlateMode('xy',0);
	setPlateMode('yz',0);
	setPlateMode('zx',0);

	var checkInput = document.getElementById('P'+plateId+'_visible');

	var mesh = canvas.getPlateMesh(plateId);

	if ( checkInput.checked ) {

		mesh.visible=true;
		mesh.geometry.elementsNeedUpdate = true;
		View_CS(0);
	}
	else {
		mesh.visible = false;
	}
}
function setPlateMode( plateId, mode )
{
	var plate = canvas.getPlate(plateId);

	var ord_XYP=0;
	var ord_YZP=0;
	var ord_ZXP=0;

	var ord_data=0;

	var Nzy =  Nz*Ny;

	var f = 0;

	var axisMax = {};
	switch( plateId ) {
		case 'xy':
			axisMax.i = Nx-1;
			axisMax.j = Ny-1;
			axisMax.ord = function( i, j ){ return (zi-1) + Nz*(j -1) + Nzy*(i -1);}
			axisMax.vf = [ [0 , Nz, Nzy ], [Nz*Ny, Nz, Nzy+Nz] ];
			break;
		case 'yz':
			axisMax.i = Ny-1;
			axisMax.j = Nz-1;
			axisMax.ord = function( i, j ){ return (j -1) + Nz*(i -1) + Nzy*(xi-1); }
			axisMax.vf = [ [ 0, 1, Nz ], [Nz, 1, Nz+1] ];
			break;
		case 'zx':
			axisMax.i = Nz-1;
			axisMax.j = Nx-1;
			axisMax.ord = function( i, j ){ return (i -1) + Nz*(yi-1) + Nzy*(j -1); }
			axisMax.vf = [ [0, Nzy, 1 ], [1, Nzy, 1+ Nzy] ];
			break;
	}

	var vertices = Vertices[0];
	
	var geometry = plate.geometry();

	var iMax = axisMax.i;
	var jMax = axisMax.j;
	switch(mode) {
		case 0:
			for (var i=1; i<=iMax; i++) {
				for (var j=1; j<=jMax; j++) {
					var ord = axisMax.ord(i, j);

					for (var k=0; k<=1; k++)
					{
						var vf = axisMax.vf[k];
						for ( var kk=0; kk <=2; kk++ ) {
							geometry.faces[f].vertexColors[kk] = new THREE.Color(
																 vertices[ord + vf[kk]].cr,
																 vertices[ord + vf[kk]].cg,
																 vertices[ord + vf[kk]].cb);
						}
						f++;
					}
				}
			}
			break;

	}
};
function movePlate( plateId ) {
	var valueInput = document.getElementById('P'+plateId+'Range');

	var V_mode = 0;

	var posInput = 'P'+plateId+'Input1';
	switch( plateId ) {
		case 'xy':
			zi = valueInput.value;

			document.getElementById(posInput).value = Vertices[0][zi-1].z;

			break;
		case 'yz':
			xi = valueInput.value;

			document.getElementById(posInput).value = Vertices[0][Nz*Ny*(xi-1)].x;

			break;
		case 'zx':
			yi = valueInput.value;

			document.getElementById(posInput).value = Vertices[0][Nz*(yi-1)].y;

			break;
	}

	setViewerShape(plateId);
	setPlateMode(plateId, V_mode);

	var plateMesh = canvas.getPlateMesh(plateId);
	plateMesh.geometry.elementsNeedUpdate = true;

	setArrows(plateId);

	View_CS(0);
}
function Sel_Draw_VR(mode)
{
	var plateId;
	var plateMesh;
	var valueInput ;
	var posInput ;

	var V_mode = 0;
	
	var i;
	
	if(mode===0 || mode===1)
	{
		
	if (document.getElementById('Pxy_visible').checked)
	{

		plateId='xy';

		valueInput = document.getElementById('P'+plateId+'Range');
		posInput = 'P'+plateId+'Input1';

		zi = valueInput.value;

		if (mode==0) { document.getElementById(posInput).value = Vertices[0][zi-1].z; }
		if (mode==1) { document.getElementById(posInput).value = Rhos[0][zi-1].z; }

		setViewerShape(plateId);
		setPlateMode(plateId, V_mode);

		plateMesh = canvas.getPlateMesh(plateId);
		plateMesh.geometry.elementsNeedUpdate = true;
	}

	if (document.getElementById('Pyz_visible').checked)
	{

		plateId='yz';

		valueInput = document.getElementById('P'+plateId+'Range');
		posInput = 'P'+plateId+'Input1';

		xi = valueInput.value;

		if (mode==0) { document.getElementById(posInput).value = Vertices[0][Nz*Ny*(xi-1)].x; }
		if (mode==1) { document.getElementById(posInput).value = Rhos[0][Nz*Ny*(xi-1)].x; }

		setViewerShape(plateId);
		setPlateMode(plateId, V_mode);

		plateMesh = canvas.getPlateMesh(plateId);
		plateMesh.geometry.elementsNeedUpdate = true;

	}

	if (document.getElementById('Pzx_visible').checked)
	{
		plateId='zx';
		valueInput = document.getElementById('P'+plateId+'Range');
		posInput = 'P'+plateId+'Input1';

		yi = valueInput.value;

		if (mode==0) { document.getElementById(posInput).value = Vertices[0][Nz*(yi-1)].y; }
		if (mode==1) { document.getElementById(posInput).value = Rhos[0][Nz*(yi-1)].y; }

		setViewerShape(plateId);
		setPlateMode(plateId, V_mode);

		plateMesh = canvas.getPlateMesh(plateId);
		plateMesh.geometry.elementsNeedUpdate = true;
	}

  }
	else
	{	
		//alert("why not");			
		
	}

	//setArrows(plateId);

}
//-----------------------------------------------
function toggleTerrain( plateId ) {

	var plate = canvas.getPlate( plateId );
	var vmodeInput = document.getElementById('P'+plateId+'_mode');
	var arrowInput = document.getElementById('P'+plateId+'_arr_visible');

	if ( vmodeInput.checked ) {
		plate.vmode( 1 );
	}
	else {
		plate.vmode( 0 );
		if ( arrowInput.checked )
			plate.visibleArrow(true);
	}

	setViewerShape( plateId );

	var mesh = plate.mesh();
	mesh.geometry.elementsNeedUpdate = true;

}

function toggleArrows( plateId ) {

	var checkInput = document.getElementById('P'+plateId+'_arr_visible');

	var plate = canvas.getPlate( plateId );
	plate.visibleArrow(checkInput.checked);
}

function show_PP_list( ai) 
{

   var objSel = document.getElementById("id_PP_list");

   objSel.options.length=0; 

   console.log("-------", ai, pseudo_type[ai][0]);

   for (i=0; i<pseudo_type[ai].length ; i++)
   {
		var objOption = document.createElement("option");       
		objOption.text = pseudo_type[ai][i];
		objOption.value = pseudo_type[ai][i];

		objSel.options.add(objOption);
  }
	
}	


function show_info(elem, Is_orb) 
{	
	if(Is_orb==1) {	F_orb =1; }
	
    var path_part = (document.getElementById('id_curr_PP_info').value).split(":")[1].trim();
	var Atom_name_temp = path_part.split(' ')[0];
	var path_part_temp = path_part.split(' ')[1];

	var data_file='/static/DB_PseudoPotentials/'+ Atom_name_temp +'/'+ path_part_temp + '/'+ Atom_name_temp +'.js3D';

	var data_string;

	const request = new XMLHttpRequest();
	const url = data_file;

	var orb_list = document.getElementById('id_orb_list');	

	if(orb_list.selectedIndex==-1) {orb_i = 1;}
	else						   {orb_i = orb_list.options[orb_list.selectedIndex].value;}


	EP_i = 0;

	request.open('GET', url, true);
	request.onload = function () 
	{  
		//data_string = request.responseText;
	};
	request.send();

	
	request.onreadystatechange = function()
	{ 
		if(request.readyState === request.DONE)
		{ 
			if(request.status === 200 || request.status === 201)
			{ 
			   console.log('성공'); 
			   data_string = request.responseText;
			   loadEPData( data_string );				  
			}  
			else
			{ console.log('실패'); 
			} 
		} 
	};

}

</script>

<script src='static/js/Plotly/plotly-2.3.0.min.js'></script>
<script>
var CS_layout = {
	title: 'Crosssection plot',
	autosize: false,
	 width: 450, height: 450,
	xaxis: {domain: [0, 1], anchor: 'y1'},
	yaxis: {domain: [0, 1], anchor: 'x1'}


  };
  
  //, tickformat : ".2f"

  function View_CS(mode)
  {	
	  if (mode ==0)
	  {

		var data_temp= V_2D;


			var data = [ {
				z: [],
				colorscale: 'Jet',
				type: 'contour',
				line:{ smoothing: 0.85 },
				xaxis: 'x1',
				yaxis: 'y1'
				}];

			data[0].z = data_temp;	

			Plotly.newPlot('CS_view', data, CS_layout);

		}

	 if (mode ==1)
	 {
		 V_2D = '';
		 var graphDiv = document.getElementById('CS_view');
		 while (graphDiv.data.length){
			Plotly.deleteTraces(graphDiv, [0]);
		}

	  }
	

  }


  var Vr_layout = {	
	autosize: false,
	margin: {	l: 20,	r: 20,	b: 20,	t: 10  },
    width: 600, height: 250,
	xaxis: { autorange: true, showgrid: true, zeroline: false, showline: true,	autotick: true,	//ticks: '',
		showticklabels: true
	},

	yaxis: { autorange: true, showgrid: true, zeroline: false, showline: true, autotick: true,	//ticks: '',
		showticklabels: true
  	}

  };




  function View_Vr(mode)
  {

	  if (mode ==0)
	  {
		

		var lq_name = ['s', 'p', 'd', 'f'];

		var trace = new Array(N_lq);
		var trace_name = new Array(N_lq);

		var data = [];

		for(var i=0; i<N_lq ;i++)
		{

			trace_name[i] = lq_name[i];  		  

			trace[i] = { x: [], y: [], name: trace_name[i], type: 'scatter' };

			trace[i].x = V1D_r[i]; trace[i].y = V1D_V[i];		  

			data.push(trace[i]);

		}

			Plotly.newPlot('Vr_view', data, Vr_layout);

	}

	if (mode ==1)
	{

		var graphDiv = document.getElementById('Vr_view');
		while (graphDiv.data.length){
			//console.log('aaa');
			Plotly.deleteTraces(graphDiv, [0]);
		}

		
	}

  }

function View_TextInfo(mode)
{
	if (mode ==0 ) {document.getElementById('id_PP_info').value = Text_Info;}
	if (mode ==1 ) {document.getElementById('id_PP_info').value = '';}
	  

}
 
</script>

	</head>
	<body style="width: 1200px; height: 1300px;">
	<h1> Pseudo potentials </h1>		
	<table style="border-collapse:collapse; border:1px gray solid;">
	<tr>
	<td colspan='1' width="1200">
	<table>"""

	
	for j in range(9):
		html_string +="<tr height='60'> "
		for i in range(18):
			html_string +="<td height='60'>"

			if pt_atom_num[j][i] !=0:

				html_string +="<button style='background-color: #" + Atom_BG_color[pt_atom_num[j][i]] + "; font-size : 17px; text-align:center; color: #ffffff ; font-weight: bold ; width:100% ; height: 30px ;margin:0;' >"
				html_string +="<sub>"  +	str(pt_atom_num[j][i]) if pt_atom_num[j][i] !=0 else ""
				html_string +="</sub>" + Atom_Name[pt_atom_num[j][i]] if pt_atom_num[j][i] !=0 else ""	
				html_string +="</button>"		 
				html_string +="<br>" if pt_atom_num[j][i] !=0 else ""
			
				html_string +="<select id='pseudo_"+Atom_Name[pt_atom_num[j][i]] + "' style='font-size:13px; color:#000000; background-color:#ffffff; width:100% ;' onchange=\"select_pseudo('" + Atom_Name[pt_atom_num[j][i]] +"','" + Atom_BG_color[pt_atom_num[j][i]] + "')\" value='test'> "
				html_string +="<option value='blank' style='background-color:#ffffff;' > </option>"					
				for p_type in pseudo_type[pt_atom_num[j][i]]:
					if p_type == '':
						html_string += "<option  selected = 'selected' value='" + p_type+"' style='color:#000000; background-color:#ffffff;'>" + p_type + "</option>"
					else : 
						html_string += "<option						value='" + p_type+"' style='color:#000000; background-color:#ffffff;'>" + p_type + "</option>"
				html_string += "</select>"
			#html_string +="</div>"
			html_string +="</td>"			
		html_string += "</tr>"
	html_string += "</table>"
	html_string += "</td>"
	html_string += "</tr>"
	html_string += "<tr>"
	html_string += "<td colspan='1'> "
	html_string += "<div id='pseu_sel' >" 
	html_string += "<input type='button' id = 'send_to_DFT' value='Send to DFT' onclick=\"send_to_DFT( )\" style='cursor:pointer; font-size:1.2em; font-weight: bold; '> "
	html_string += "</div> "
	html_string +="</td>"
	html_string +="</tr>"
	html_string +="<tr><td>"
	html_string +="    <table border='1'><tr>"
	html_string +="           <td style=' height: 250px; width: 400px;'>" 
	html_string +="           <input id='id_curr_PP_info' type='text' style=' height: 30px; width: 400px; font-weight: bold; font-size:15px'  autocomplete='off'>"
	html_string +="           <textarea id='id_PP_info' style=' height: 220px; width: 400px; font-weight: bold; font-size:15px'; >	</textarea> "
	html_string +="           </td>"
	html_string +="           <td style=' height: 250px; width: 800px;'> <div id='Vr_view' style='width:100%;margin: auto;'></div> </td>"
	html_string +="    </tr></table>"

	html_string +="</td></tr>"
	
	html_string +="<tr><td>"
	html_string +="    <table border='1'><tr>"

	html_string +="           <td style=' width: 35%;'> <div id='CS_view' style='width:100%;margin: auto;'></div> </td>"
	html_string +="           <td style=' width: 40%;'> <div id='container' style='width:100%;margin: auto;'> </div> </td>"
	html_string +="           <td style=' width: 25%;'> <div style='width:100%;margin: auto;'> Orbitals "
	html_string +="           <select id = 'id_orb_list' style='font-size:11px; color:#000000; background-color:#ffffff;' onchange='show_info(0,1)'>"
	html_string +="           </select>"	
	html_string +="           </div>"	

	html_string +="           <table> "
	html_string +="           <tr> "
	html_string +="           <td> Value level </td> "
	html_string +="           <td> "
	html_string +="           <input id='MC_Range' type='range' style='width:90%' value='0.5' min='0' max='1' step='0.1' step='0.1' oninput=\"MC_effect()\" onmouseover=\"OFF_MouseMove()\" onmouseleave=\"ON_MouseMove()\" autocomplete='off'/>"
	html_string +="           </td> "
	html_string +="           <td> <input id='MC_Range_text' type='text' value='0.5' size='1' height='1' autocomplete='off'> </td> "

	html_string +="           </tr> "
	html_string +="           </table> "

	html_string +="           <table> "
	html_string +="           <tr bgcolor='#eeeeee' > "
	html_string +="           <td> <SPAN style='font-size: 10pt'> Viewer </SPAN> </td>"
	html_string +="           <td> <SPAN style='font-size: 10pt'> density </SPAN> </td>"
	html_string +="           <td> <SPAN style='font-size: 10pt'> terrain </SPAN> </td>"
	html_string +="           <td> <SPAN style='font-size: 10pt'> E field </SPAN> </td> "
	html_string +="            </tr>"
	html_string +="           <tr> "
	html_string +="           <td> <SPAN style='font-size: 10pt'> Pxy </SPAN> </td>"
	html_string +="           <td> <input id='Pxy_visible' type='checkbox' onclick=\"toggleView( 'xy' )\"> </td>"
	html_string +="           <td> <input id='Pxy_mode' type='checkbox' onclick=\"toggleTerrain('xy')\"> </td> "
	html_string +="           <td> <input id='Pxy_arr_visible' type='checkbox' onclick=\"toggleArrows('xy')\"> </td> "
	html_string +="           </tr>"

	html_string +="           <tr>"
	html_string +="           <td> <SPAN style='font-size: 10pt'> Pyz </SPAN> </td>"
	html_string +="           <td> <input id='Pyz_visible' type='checkbox' onclick=\"toggleView( 'yz' )\"> </td>"
	html_string +="           <td> <input id='Pyz_mode' type='checkbox' onclick=\"toggleTerrain('yz')\"> </td> "
	html_string +="           <td> <input id='Pyz_arr_visible' type='checkbox' onclick=\"toggleArrows('yz')\"> </td> "
	html_string +="            </tr>"
	html_string +="           <tr> "
	html_string +="           <td> <SPAN style='font-size: 10pt'> Pzx </SPAN> </td>"
	html_string +="           <td> <input id='Pzx_visible' type='checkbox' onclick=\"toggleView( 'zx' )\"> </td> "
	html_string +="           <td> <input id='Pzx_mode' type='checkbox' onclick=\"toggleTerrain('zx')\"> </td> "
	html_string +="           <td> <input id='Pzx_arr_visible' type='checkbox' onclick=\"toggleArrows('zx')\"> </td> "
	html_string +="           </tr> "
	html_string +="           </table> "
	html_string +="           <table style=\"width:95%\"> "
	html_string +="           <tr> "
	html_string +="           <td colspan='1'> <SPAN style='font-size: 10pt'> Viewer </SPAN> </td>"
	html_string +="           <td colspan='1' bgcolor='#eeeeee'> <SPAN style='font-size: 10pt'> position </SPAN></td>"
	html_string +="           <td bgcolor='#eeeeee'> </td>"
	html_string +="           </tr> "
	html_string +="           <tr> "
	html_string +="           <td colspan='1'> <SPAN style='font-size: 10pt'> Pxy </SPAN> </td>"
	html_string +="           <td colspan='1'> <input id='PxyRange' type='range' style='width:90%' value='0' min='1' max='1' step='1' step='1' oninput=\"movePlate('xy')\" onmouseover=\"OFF_MouseMove()\" onmouseleave=\"ON_MouseMove()\" autocomplete='off'/> </td>"
	html_string +="           <td> <input id='PxyInput1' type='text' value='0' size='3' height='1' autocomplete='off'> </td>"
	html_string +="           </tr> "
	html_string +="           <tr> "
	html_string +="           <td colspan='1'> <SPAN style='font-size: 10pt'> Pyz </SPAN> </td>"
	html_string +="           <td colspan='1'> <input id='PyzRange' type='range' style='width:90%' value='0' min='1' max='1' step='1' oninput=\"movePlate('yz')\" onmouseover=\"OFF_MouseMove()\" onmouseleave=\"ON_MouseMove()\" autocomplete='off'/> </td>"
	html_string +="           <td> <input id='PyzInput1' type='text' value='0' size='3' height='1' autocomplete='off'> </td>"
	html_string +="           </tr>"
	html_string +="           <tr> "
	html_string +="           <td colspan='1'> <SPAN style='font-size: 10pt'> Pzx </SPAN> </td>"
	html_string +="           <td colspan='1'> <input id='PzxRange' type='range' style='width:90%' value='0' min='1'max='1' step='1' oninput=\"movePlate('zx')\" onmouseover=\"OFF_MouseMove()\" onmouseleave=\"ON_MouseMove()\" autocomplete='off'/> </td>"
	html_string +="           <td> <input id='PzxInput1' type='text' value='0' size='3' height='1' autocomplete='off'> </td>"
	html_string +="           </tr>"
	html_string +="           </table>"

	html_string += "	</td>"
	html_string += "    </tr></table>"

	html_string += "</td></tr>"
	
	html_string += "</table>"
	html_string += "</body>"
	html_string += "</html>"
	
	
	#command = "rm -f " + "static/htmls/pseudo_ini_*.html"; os.system(command)
	#PP_url = "static/htmls/" + "pseudo_ini_"+str(int(time.time()*10000000))+".html"
	PP_url = "static/htmls/" + "Pseudo_ini.html"

	p_file = open(PP_url, 'w')
	p_file.write(html_string)
	p_file.close()

	return PP_url

def create_NewLocalDB_ini():

	# load user names ----------------------------------

	client = MongoClient('localhost', 27017)

	db = client['NewDBmat']	
	
	#Users = db.client["DBmat"]["Group"].distinct("User")


	#Using pymongo and connection named 'db' reduce( lambda all_keys, rec_keys: all_keys | set(rec_keys), map(lambda d: d.keys(), db.things.find()), set() )
	
	# read orbs------------------------------------------
		

	command = "cp static/NewLocalDB_ini0.html templates/NewLocalDB_ini.html"

	os.system(command)			

	return 

def create_LocalDB():

	# load user names ----------------------------------

	client = MongoClient('localhost', 27017)

	db = client['DBmat']	
	
	Users = db.client["DBmat"]["Group"].distinct("User")

 
	# read ------------------------------------------
	Read_html = open('static/LocalDB_ini.html', 'r')
	
	save_string =''
	
	while True:
		html_line = Read_html.readline()  
		if not html_line : break
		
		change_string = "<body onload ='Initialize();'"
		if change_string in html_line: html_line = html_line.replace(change_string, "<body onload ='Initialize({{selected_i}});'")

		
		change_string = '<option id="---op1---"> </option>'

		if change_string in html_line: 
			replace_string  = '{% for key_value in field_key[1] %}'
			replace_string += '<option {% if (loop.index0==selected_i[1]|int) %} selected {% endif %} > {{key_value}} </option>'
			replace_string += '{% endfor %}'
			html_line = html_line.replace(change_string, replace_string)

		change_string = '<option id="---op2---"> </option>'

		if change_string in html_line: 
			replace_string  = '{% for key_value in field_key[2] %}'
			replace_string += '<option {% if (loop.index0==selected_i[2]|int) %} selected {% endif %} > {{key_value}} </option>'
			replace_string += '{% endfor %}'
			html_line = html_line.replace(change_string, replace_string)

		change_string = '<option id="---op3---"> </option>'

		if change_string in html_line: 
			replace_string  = '{% for key_value in field_key[3] %}'
			replace_string += '<option {% if (loop.index0==selected_i[3]|int) %} selected {% endif %} > {{key_value}} </option>'
			replace_string += '{% endfor %}'
			html_line = html_line.replace(change_string, replace_string)

		change_string = '<option id="---op4---"> </option>'

		if change_string in html_line: 
			replace_string  = '{% for key_value in field_key[4] %}'
			replace_string += '<option {% if (loop.index0==selected_i[4]|int) %} selected {% endif %} > {{key_value}} </option>'
			replace_string += '{% endfor %}'
			html_line = html_line.replace(change_string, replace_string)

		'''
		change_string = '<option id="---op5---"> </option>'

		if change_string in html_line: 
			replace_string  = '{% for key_value in field_key[5] %}'
			replace_string += '<option {% if (loop.index0==selected_i[5]|int) %} selected {% endif %} > {{key_value}} </option>'
			replace_string += '{% endfor %}'
			html_line = html_line.replace(change_string, replace_string)
		'''
		#selected_i
		#var sel = document.getElementById("exam"); for(var i=0; i<sel.length; i++){ if(sel[i].value==2){ sel[i].selected = true; }
			
		#field_key

		change_string = "name = 'Info_path' value=''"
		if change_string in html_line: html_line = html_line.replace(change_string, "name = 'Info_path' value='{{DB_key['Info_path']}}'")

		change_string = "name = 'search_text' value=''"
		if change_string in html_line: html_line = html_line.replace(change_string, "name = 'search_text' value='{{DB_key['search_text']}}'")


		change_string = "id='id_iframe2'"
		if change_string in html_line: html_line = html_line.replace(change_string, """src="{{ url_for('static', filename=center_url) }}" """)

		change_string = "items has been found"
		if change_string in html_line: html_line = html_line.replace(change_string, """ {{N_db}} items has been found""")

		

		save_string += html_line 
				  
	Read_html.close()

	  
	html_url = "templates/LocalDB.html"

	p_file = open(html_url, 'w')
	p_file.write(save_string)
	p_file.close()

	return html_url

#--------------------------------------------------------------------------
def create_DFT():	
		
	command = "cp static/Run_DFT_ini.html templates/Run_DFT.html"

	os.system(command)			

	return 0


#--------------------------------------------------------------------------
def create_Run_button_html():

	print("create_Run_button_html")

	command = "cp static/Run_button.html templates/Run_button.html"

	os.system(command)			

	return 0
#--------------------------------------------------------------------------
def create_NanoCore_button_html():

	command = "cp static/NanoCore_button.html templates/NanoCore_button.html"

	os.system(command)			

	return 0

#-------------------------------------------------------------------------
def create_Fit_button_html():

	command = "cp static/Fit_button.html templates/Fit_button.html"

	os.system(command)			

	return 0
#--------------------------------------------------------------------------
def create_Login_SIESTA():

	# read ------------------------------------------
	Read_html = open('static/Login_SIESTA_ini.html', 'r')

	save_string =''

	while True:
		html_line = Read_html.readline()  
		if not html_line : break
		
		change_string = "<body>"
		if change_string in html_line: html_line = html_line.replace(change_string, "<body onload='Initialize({{F_initial}})'")

		change_string = "id='id_id_text2'"
		if change_string in html_line: html_line = html_line.replace(change_string, "id='id_id_text2' value='{{id_text}}' ")

		change_string = "id='id_pass_text2'"
		if change_string in html_line: html_line = html_line.replace(change_string, "id='id_pass_text2' value='{{pass_text}}' ")
		
		save_string += html_line
	
	Read_html.close()

	print(save_string)

	
	html_url = "templates/Login_SIESTA.html"

	p_file = open(html_url, 'w')
	p_file.write(save_string)
	p_file.close()
	

	return 0	
#-------------------------------------------------------------------------
def create_Dropzone_html():

	command = "cp static/Dropzone0.html templates/Dropzone.html"

	os.system(command)			

	return 0

#-------------------------------------------------------------------------
def create_Login_DB_html():

	command = "cp static/Login_DB_ini.html templates/Login_DB.html"

	os.system(command)			

	return 0	

#-------------------------------------------------------------------------
def create_Logout_DB_html():

	command = "cp static/Logout_DB_ini.html templates/Logout_DB.html"

	os.system(command)			

	return 0		
#-------------------------------------------------------------------------
def create_Register_DB_html():

	command = "cp static/Register_DB_ini.html templates/Register_DB.html"

	os.system(command)			

	return 0			
#-------------------------------------------------------------------------
def create_DB_search():

	command = "cp static/DB_search.html templates/DB_search.html"

	os.system(command)			

	return 0				
#-------------------------------------------------------------------------
def create_DB_updateWeb():

	command = "cp static/DB_updateWeb.html templates/DB_updateWeb.html"

	os.system(command)			

	return 0					