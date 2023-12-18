//	Program developed by
//	
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut
//
//
//

#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include "cargar-triangulo.h"
#include <math.h>
#include "obj.h"

#define ANGELUA 0.2
#define DESPL 5
#define SCALE 1.1

/*
typedef struct mlist
    {
    double m[16];
    struct mlist *hptr;
    } mlist;
*/
typedef struct triobj
    {
    hiruki *triptr;
    int num_triangles;
    mlist *mptr;
    struct triobj *hptr;
    double modelview[16];
    int is_cam;
    unsigned char * rgb;
    int has_color;
    } triobj;

typedef struct camera_obj
{
    double m_obj[16];
    double m_esa[16];
    struct camera_obj *prev_camera;
} camera_obj;
// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy;

int indexx;
hiruki *triangulosptr;
object3d *foptr;
object3d *sel_ptr;
triobj *cam_obj_ptr;
camera_obj *camera;
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;
int resize;
int perspective;
int back_culling;
int show_normal_v;
int analisis_mode;
int camera_selected; //0= Object selected, 1= Camera selected
int object_as_camera; //0= Camera is not an Object, 1= Camera is an object
double modelview[16];
double modelview2[16];
double camera_by_mx[16];
double projection_mx[16];

char fitxiz[100];

void update_cam_obj(){
    cam_obj_ptr->mptr->m[3]=camera->m_obj[3];
    cam_obj_ptr->mptr->m[7]=camera->m_obj[7];
    cam_obj_ptr->mptr->m[11]=camera->m_obj[11];
}

void matrix_calc(double * result, double * left_mx, double * right_mx){
    double bat;
    int i,j,k;
    for(i=0; i<4;i++){
        for(j=0;j<4;j++){
            bat=0.0;
            for(k=0; k<4;k++){
                bat+=left_mx[i*4+k]*right_mx[k*4+j];
            }
            result[i*4+j] = bat;
        }
    }

}

void mxp(punto *pptr, double m[16], punto p)
{
//pptr->x = p.x;
pptr->x = (m[0]*p.x)+(m[1]*p.y)+(m[2]*p.z)+m[3];
pptr->y = (m[4]*p.x)+(m[5]*p.y)+(m[6]*p.z)+m[7];
pptr->z = (m[8]*p.x)+(m[9]*p.y)+(m[10]*p.z)+m[11];
pptr->u = p.u;
pptr->v = p.v;
}

int mxprojection(punto * point, double mx[16], punto og_pnt){
    double nx,ny,nz,w;

    nx = (mx[0]*og_pnt.x)+(mx[1]*og_pnt.y)+(mx[2]*og_pnt.z)+mx[3];
    ny = (mx[4]*og_pnt.x)+(mx[5]*og_pnt.y)+(mx[6]*og_pnt.z)+mx[7];
    nz = (mx[8]*og_pnt.x)+(mx[9]*og_pnt.y)+(mx[10]*og_pnt.z)+mx[11];
    w = (mx[12]*og_pnt.x)+(mx[13]*og_pnt.y)+(mx[14]*og_pnt.z)+mx[15];
    //printf("%f\n",w);

    if (w<=0) return 0;

    nx = nx*500/w;
    ny = ny*500/w;
    nz = nz*-500/w;

    point->x = nx;
    point->y = ny;
    point->z = nz;

    return 1;
}

double calc_normal(double * vector){
    double answ =  sqrt((vector[0]*vector[0])+(vector[1]*vector[1])+(vector[2]*vector[2]));
    if (answ<0){
        answ = -answ;
    }
    return answ;
}

void normalize(double * result, double * vector, double normal){
    result[0] = vector[0]/normal;
    result[1] = vector[1]/normal;
    if (result[1]<0) result[1] = -result[1];
    result[2] = vector[2]/normal;
    if (result[2]<0) result[2] = -result[2];
}

void m_esa_calc(camera_obj *cam){
    //CHANGES DONE
    double e_x,e_y,e_z,x_c,y_c,z_c,e;
    e_x = (cam->m_obj[3]);
    e_y = (cam->m_obj[7]);
    e_z = (cam->m_obj[11]);
    x_c = cam->m_obj[0]*(-e_x) + cam->m_obj[4]*(-e_y) + cam->m_obj[8]*(-e_z);
    y_c = cam->m_obj[1]*(-e_x) + cam->m_obj[5]*(-e_y) + cam->m_obj[9]*(-e_z);
    z_c = cam->m_obj[2]*(-e_x) + cam->m_obj[6]*(-e_y) + cam->m_obj[10]*(-e_z);
    cam->m_esa[0]= cam->m_obj[0]; cam->m_esa[1]=cam->m_obj[4]; cam->m_esa[2]=cam->m_obj[8]; cam->m_esa[3]=x_c;
    cam->m_esa[4]= cam->m_obj[1]; cam->m_esa[5]=cam->m_obj[5]; cam->m_esa[6]=cam->m_obj[9]; cam->m_esa[7]=y_c;
    cam->m_esa[8]=cam->m_obj[2]; cam->m_esa[9]=cam->m_obj[6]; cam->m_esa[10]=cam->m_obj[10]; cam->m_esa[11]=z_c;
    cam->m_esa[12]=cam->m_obj[3]; cam->m_esa[13]=cam->m_obj[7]; cam->m_esa[14]=cam->m_obj[11]; cam->m_esa[15]=cam->m_obj[15];
}

void calc_v_normal(hiruki * tri, punto p1, punto p2, punto p3){
    float df1[3];
    float df2[3];
    double norm_len;
    df1[0] = p2.x - p1.x; df1[1] = p2.y - p1.y; df1[2] = p2.z - p1.z;
    df2[0] = p3.x - p1.x; df2[1] = p3.y - p1.y; df2[2] = p3.z - p1.z;
    tri->bek_normal.x = (df1[1]*df2[2])-(df1[2]*df2[1]);
    tri->bek_normal.y = (df1[2]*df2[0])-(df1[0]*df2[2]);
    tri->bek_normal.z = (df1[0]*df2[1])-(df1[1]*df2[0]);
    norm_len = sqrt(pow(tri->bek_normal.x,2)+pow(tri->bek_normal.y,2)+pow(tri->bek_normal.z,2));
    tri->bek_normal.x /= norm_len;
    tri->bek_normal.y /= norm_len;
    tri->bek_normal.z /= norm_len;
}

int is_facing(hiruki * tri, punto p){
    float cam_ray[3];
    float tri_mod,cam_mod,mod;
    double vn,ang;
    punto p_cam;

    p_cam.x = camera->m_obj[3];p_cam.y=camera->m_obj[7];p_cam.z=camera->m_obj[11];
    //mxp(&p_cam,modelview,p_cam);
    
    // Show Point to Camera Line
    /*glBegin(GL_LINES);
    glVertex3d(p.x,p.y,p.z); //tri->p1.x,tri->p1.y,tri->p1.z
    glVertex3d(p_cam.x,p_cam.y,p_cam.z);
    glEnd();*/
    
    

    cam_ray[0]= p_cam.x-p.x; //camera->m_obj[3]-p.x
    cam_ray[1]= p_cam.y-p.y;
    cam_ray[2]= p_cam.z-p.z;
    cam_mod = sqrt(pow(cam_ray[0],2)+pow(cam_ray[1],2)+pow(cam_ray[2],2));
    vn = (cam_ray[0]*tri->bek_normal.x) + (cam_ray[1]*tri->bek_normal.y) + (cam_ray[2]*tri->bek_normal.z);
    
    tri_mod = sqrt(pow(tri->bek_normal.x,2)+pow(tri->bek_normal.y,2)+pow(tri->bek_normal.z,2));
    cam_mod = sqrt(pow(cam_ray[0],2)+pow(cam_ray[1],2)+pow(cam_ray[2],2));
    mod = tri_mod * cam_mod;
    ang = vn/mod;

    if (ang<0){
        return 0;
    }
    return 1;
}

void vertex_2_point(punto * p, vertex v){
    p->u = v.u;
    p->v = v.v;
    p->x = v.coord.x;
    p->y = v.coord.y;
    p->z = v.coord.z;
}

void objektuari_aldaketa_sartu_ezk(double m[16])
{
}
void objektuari_aldaketa_sartu_esk(double m[16])
{
}

unsigned char * color_textura(float u, float v)
{
int desplazamendua,ix_x,ix_y;
char * lag;

if(u<0) u=0;
if(u>1) u=1;
if(v<0) v=0;
if(v>1) v=1;

ix_x = (int)(u*(dimx-1));
ix_y = (int)((1-v)*(dimy-1));

desplazamendua = ix_y*dimx+ix_x;

lag = (unsigned char *)bufferra;  // pixel on the left and top
return(lag+3*desplazamendua);
}

void  dibujar_linea_z(int linea,float c1x, float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v)
{
float xkoord,zkoord;
float u,v,t,dt,cx,cz,cu,cv;
unsigned char r,g,b;
unsigned char *colorv;

//printf("%f",c1u);

glBegin( GL_POINTS );

if (c2x == c1x){
    dt = 1;
}else{
    dt = 1/(c2x-c1x);
}

for (xkoord = c1x, zkoord=c1z, u=c1u, v=c1v, t=1; xkoord<=c2x; xkoord++, t= t-dt)
    {
        colorv = color_textura(u,v);
    }
    r= colorv[0];
    g=colorv[1];
    b=colorv[2];  
    //r = 0; g=0; b=0;


    glColor3ub(r,g,b);
    glVertex3f(xkoord, linea, zkoord );

    u = t*c1u + (1-t)*c2u;
    v = t*c1v + (1-t)*c2v;
    zkoord = t*c1z + (1-t)*c2z;
    //zkoord = 0.0;
glEnd();
}

void print_matrizea(double *m)
{
int i;

for (i = 0;i<4;i++)
   printf("%lf, %lf, %lf, %lf\n",m[i*4],m[i*4+1],m[i*4+2],
                                 m[i*4+3]);
}


// TODO
// aurrerago egitekoa
// para más adelante

void ebakidura_kalk(punto *gptr, punto *bptr, int h, punto *emptr){

//emptr->y kalkulatu
emptr->y=(float)h;

float dify2 = h-bptr->y;
float dify = gptr->y - bptr-> y;
if (dify >0){
    //emptr->x kalkulatu
    float difx = gptr->x-bptr->x;
    emptr->x= ((dify2/dify)*difx)+bptr->x;

    //emptr->z kalkulatu
    float difz = gptr->z-bptr->z;
    emptr->z= ((dify2/dify)*difz)+bptr->z;

    //emptr->u kalkulatu
    float difu = gptr->u-bptr->u;
    emptr->u= ((dify2/dify)*difu)+bptr->u;

    //emptr->v kalkulatu
    float difv = gptr->v-bptr->v;
    emptr->v= ((dify2/dify)*difv)+bptr->v;
}else{
    emptr->x = bptr->x;
    emptr->z = bptr->z;
    emptr->u = bptr->u;
    emptr->v = bptr->v;
}

}

//THE OLD
/*
void dibujar_triangulo(triobj *optr, int i)
{
hiruki *tptr;
punto *pgoiptr, *pbeheptr, *perdiptr, *middle;
float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v,tald,qald,sald;
int linea,t,q,s,xpersp,ypersp,zpersp;
float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
punto p1,p2,p3,p_sar,p_irt,p_help,pp1,pp2,pp3;

if (i >= optr->num_triangles) return;
tptr = optr->triptr+i;


matrix_calc(modelview,camera->m_esa,optr->mptr->m);
mxp(&p1,modelview,tptr->p1);
mxp(&p2,modelview,tptr->p2);
mxp(&p3,modelview,tptr->p3); //optr->mptr->m //optr->modelview


if (perspective){
    //matrix_calc(modelview2,projection_mx,modelview);
    xpersp=mxprojection(&p1,projection_mx,p1);
    ypersp=mxprojection(&p2,projection_mx,p2);
    zpersp=mxprojection(&p3,projection_mx,p3);
    if (xpersp==0 || ypersp==0 || zpersp == 0){
        return;
    }
}

calc_v_normal(tptr, p1,p2,p3);
//is_facing(tptr);

if (show_normal_v) {
    glBegin(GL_LINES);
    glVertex3d(p1.x,p1.y,p1.z);
    glVertex3d(p1.x+(100*tptr->bek_normal.x),p1.y+(100*tptr->bek_normal.y),p1.z+(100*tptr->bek_normal.z));
    glEnd();
}

int is_facing_val = is_facing(tptr, p1);

if (lineak == 1){
    if (is_facing_val) {
        glBegin(GL_POLYGON);
        glColor3ub(255,255,255);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }else if (!is_facing_val && back_culling){
        glBegin(GL_POLYGON);
        glColor3ub(255,0,0);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }
    return;
}
        
        //  else 
  
if(p1.y >=p2.y) {pgoiptr = &(p1); pbeheptr = &(p2);}
else{pgoiptr = &(p2); pbeheptr = &(p1);}
if(p3.y >= pgoiptr->y) {perdiptr = pgoiptr;pgoiptr = &(p3);}
else if(p3.y<pbeheptr->y){perdiptr=pbeheptr; pbeheptr = &(p3);}
else{perdiptr=&(p3);};

if(perdiptr->y == pbeheptr->y && perdiptr->x>=pbeheptr->x){
    middle = perdiptr;
    perdiptr = pbeheptr;
    pbeheptr = middle;
}
if(perdiptr->y==pgoiptr->y && pgoiptr->x<=perdiptr->x){
    middle = perdiptr;
    perdiptr = pgoiptr;
    pgoiptr = middle;

}

for(int i= pgoiptr->y;i>perdiptr->y;i--){
    ebakidura_kalk(pgoiptr, perdiptr, i, &p_sar);
    ebakidura_kalk(pgoiptr, pbeheptr, i, &p_irt);
    if(p_sar.x>p_irt.x){p_help = p_sar;p_sar = p_irt;p_irt = p_help;}
    dibujar_linea_z(i, p_sar.x, p_sar.z, p_sar.u, p_sar.v, p_irt.x, p_irt.z, p_irt.u, p_irt.v,optr->has_color,optr->rgb);
}
for(int i=perdiptr->y;i>pbeheptr->y;i--){
    ebakidura_kalk(perdiptr, pbeheptr, i, &p_sar);
    ebakidura_kalk(pgoiptr, pbeheptr, i, &p_irt);
    if(p_sar.x>p_irt.x){p_help = p_sar;p_sar = p_irt;p_irt = p_help;}
    dibujar_linea_z(i, p_sar.x, p_sar.z, p_sar.u, p_sar.v, p_irt.x, p_irt.z, p_irt.u, p_irt.v,optr->has_color,optr->rgb);
}
}
*/

//THE NEW
void dibujar_triangulo(vertex vx1, vertex vx2, vertex vx3, object3d * obj)
{
hiruki *tptr;
punto *pgoiptr, *pbeheptr, *perdiptr, *middle;
float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v,tald,qald,sald;
int linea,t,q,s,xpersp,ypersp,zpersp;
float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
punto p1,p2,p3,p_sar,p_irt,p_help;
vertex_2_point(&p1,vx1);
vertex_2_point(&p2,vx2);
vertex_2_point(&p3,vx3);

matrix_calc(modelview,camera->m_esa,obj->mptr->m);
mxp(&p1,modelview,p1);
mxp(&p2,modelview,p2);
mxp(&p3,modelview,p3); //optr->mptr->m //optr->modelview


if (perspective){
    //matrix_calc(modelview2,projection_mx,modelview);
    xpersp=mxprojection(&p1,projection_mx,p1);
    ypersp=mxprojection(&p2,projection_mx,p2);
    zpersp=mxprojection(&p3,projection_mx,p3);
    if (xpersp==0 || ypersp==0 || zpersp == 0){
        return;
    }
}

//calc_v_normal(tptr, p1,p2,p3);
//is_facing(tptr);
/*
if (show_normal_v) {
    glBegin(GL_LINES);
    glVertex3d(p1.x,p1.y,p1.z);
    glVertex3d(p1.x+(100*tptr->bek_normal.x),p1.y+(100*tptr->bek_normal.y),p1.z+(100*tptr->bek_normal.z));
    glEnd();
}

int is_facing_val = is_facing(tptr, p1);
*/
int is_facing_val = 1;
if (lineak == 1){
    if (is_facing_val) {
        glBegin(GL_POLYGON);
        glColor3ub(255,255,255);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }else if (!is_facing_val && back_culling){
        glBegin(GL_POLYGON);
        glColor3ub(255,0,0);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }
    return;
}
       
        //  else 
  
if(p1.y >=p2.y) {pgoiptr = &(p1); pbeheptr = &(p2);}
else{pgoiptr = &(p2); pbeheptr = &(p1);}
if(p3.y >= pgoiptr->y) {perdiptr = pgoiptr;pgoiptr = &(p3);}
else if(p3.y<pbeheptr->y){perdiptr=pbeheptr; pbeheptr = &(p3);}
else{perdiptr=&(p3);};

if(perdiptr->y == pbeheptr->y && perdiptr->x>=pbeheptr->x){
    middle = perdiptr;
    perdiptr = pbeheptr;
    pbeheptr = middle;
}
if(perdiptr->y==pgoiptr->y && pgoiptr->x<=perdiptr->x){
    middle = perdiptr;
    perdiptr = pgoiptr;
    pgoiptr = middle;

}

for(int i= pgoiptr->y;i>perdiptr->y;i--){
    ebakidura_kalk(pgoiptr, perdiptr, i, &p_sar);
    ebakidura_kalk(pgoiptr, pbeheptr, i, &p_irt);
    if(p_sar.x>p_irt.x){p_help = p_sar;p_sar = p_irt;p_irt = p_help;}
    dibujar_linea_z(i, p_sar.x, p_sar.z, p_sar.u, p_sar.v, p_irt.x, p_irt.z, p_irt.u, p_irt.v);
}
for(int i=perdiptr->y;i>pbeheptr->y;i--){
    ebakidura_kalk(perdiptr, pbeheptr, i, &p_sar);
    ebakidura_kalk(pgoiptr, pbeheptr, i, &p_irt);
    if(p_sar.x>p_irt.x){p_help = p_sar;p_sar = p_irt;p_irt = p_help;}
    dibujar_linea_z(i, p_sar.x, p_sar.z, p_sar.u, p_sar.v, p_irt.x, p_irt.z, p_irt.u, p_irt.v);
}
}

//THE OLD 
/*
static void marraztu(void)
{
float u,v;
int i,j;
triobj *auxptr;

unsigned char* colorv;
unsigned char r,g,b;


  // marrazteko objektuak behar dira
  // no se puede dibujar sin objetos
if (foptr ==0) return;

// clear viewport...
if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
    else 
      {
      if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
      }

glMatrixMode(GL_PROJECTION);
glLoadIdentity();
glOrtho(-500.0, 500.0, -500.0, 500.0, -500, 500.0);


//triangulosptr = sel_ptr->triptr;
if (objektuak == 1)
    {
    if (denak == 1)
        {
        for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
            for (i =0; i < auxptr->num_triangles; i++)
                {
                dibujar_triangulo(auxptr,i);
                }
            }
        }
      else
        {
        for (i =0; i < sel_ptr->num_triangles; i++)
            {
            dibujar_triangulo(sel_ptr,i);
            }
        }
    }
  else
    {
     dibujar_triangulo(sel_ptr,indexx);
    }
glFlush();
}
*/

//THE NEW
static void marraztu(void){
    object3d * next_object;
    face next_face;
    vertex vertex1,vertex2,vertex3;
    int i,j;
    if (foptr ==0) return;

    // clear viewport...
    if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
        else 
        {
        if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
        }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-500.0, 500.0, -500.0, 500.0, -500, 500.0);

    if (objektuak){
        if (denak){
            for(next_object = foptr; next_object!=0; next_object=next_object->hptr){
                for (i = 0; i<=next_object->num_faces-1;i++){
                    next_face = next_object->face_table[i];
                    vertex1 = next_object->vertex_table[next_face.vertex_ind_table[0]];
                    for (j=2; j<=next_face.num_vertices-1;j++){
                        vertex2 = next_object->vertex_table[next_face.vertex_ind_table[j-1]];
                        vertex3 = next_object->vertex_table[next_face.vertex_ind_table[j]];
                        
                        dibujar_triangulo(vertex1,vertex2,vertex3,next_object);
                        //TODO: DRAW TRIANGLE
                    }
                }
            }
        }
    }
}

void read_from_file(char *fitx)
{
int i,retval;
object3d *optr;

    //printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (object3d *)malloc(sizeof(object3d));
    //retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &(optr->rgb));
    retval = read_wavefront(fitx,optr);
    if (retval !=0) 
         {
         printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n",fitxiz);
         free(optr);
         }
       else
         {
            /*
         triangulosptr = optr->triptr;
         //printf("objektuaren matrizea...\n");
         optr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         optr->mptr->m[0] = 1.0;
         optr->mptr->m[5] = 1.0;
         optr->mptr->m[10] = 1.0;
         optr->mptr->m[15] = 1.0;
         optr->mptr->hptr = 0;
         optr->is_cam=0;*/
         //printf("objektu zerrendara doa informazioa...\n");
         optr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         optr->mptr->m[0] = 1.0;
         optr->mptr->m[5] = 1.0;
         optr->mptr->m[10] = 1.0;
         optr->mptr->m[15] = 1.0;
         optr->mptr->hptr = 0;
         optr->hptr = foptr;
         foptr = optr;
         sel_ptr = optr;
         }
     printf("datuak irakurrita\nLecura finalizada\n");
}

void convert_object_to_camera(){
    camera_obj* new_cam = (camera_obj *)malloc(sizeof(camera_obj));
    new_cam->m_obj[0] = sel_ptr->mptr->m[0]; new_cam->m_obj[1] = sel_ptr->mptr->m[1]; new_cam->m_obj[2] = sel_ptr->mptr->m[2]; new_cam->m_obj[3] = sel_ptr->mptr->m[3];
    new_cam->m_obj[4] = sel_ptr->mptr->m[4]; new_cam->m_obj[5] = sel_ptr->mptr->m[5]; new_cam->m_obj[6] = sel_ptr->mptr->m[6]; new_cam->m_obj[7] = sel_ptr->mptr->m[7];
    new_cam->m_obj[8] = sel_ptr->mptr->m[8]; new_cam->m_obj[9] = sel_ptr->mptr->m[9]; new_cam->m_obj[10] = sel_ptr->mptr->m[10]; new_cam->m_obj[11] = sel_ptr->mptr->m[11];
    new_cam->m_obj[12] = sel_ptr->mptr->m[12]; new_cam->m_obj[13] = sel_ptr->mptr->m[13]; new_cam->m_obj[14] = sel_ptr->mptr->m[14]; new_cam->m_obj[15] = sel_ptr->mptr->m[15];
    m_esa_calc(new_cam);
    new_cam->prev_camera = camera;
    camera = new_cam;
    sel_ptr->is_cam = 1;
}

void update_camera(double * mx){
    camera_obj* new_cam = (camera_obj *)malloc(sizeof(camera_obj));
    double * cam = camera->m_obj; //Current Camera Matrix
    matrix_calc(&(new_cam->m_obj[0]), &(cam[0]) ,&(mx[0]));
    m_esa_calc(new_cam);
    new_cam->prev_camera = camera;
    camera = new_cam;
}

void form_matrix_x(double *mx,int dir){
    switch (aldaketa)
    {
    case 'r':
        mx[0]=1; mx[1]= 0; mx[2]=0; mx[3]= 0;
        mx[4]=0; mx[5]= cos(dir*ANGELUA); mx[6]= -sin(dir*ANGELUA); mx[7]= 0;
        mx[8]=0; mx[9]= sin(dir*ANGELUA); mx[10]= cos(dir*ANGELUA); mx[11]= 0;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
        break;
    
    case 't':
        mx[0]=1; mx[1]=0; mx[2]=0; mx[3]=dir*DESPL;
        mx[4]=0; mx[5]=1; mx[6]=0; mx[7]=0;
        mx[8]=0; mx[9]=0; mx[10]=1; mx[11]=0;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
        break;

    default:
        break;
    }
}

void form_matrix_y(double *mx, int dir){
    switch (aldaketa)
    {
    case 'r':
        mx[0]=cos(dir*ANGELUA); mx[1]=0; mx[2] =sin(dir*ANGELUA); mx[3]=0;
        mx[4]=0; mx[5]=1; mx[6]=0; mx[7]=0;
        mx[8]=-sin(dir*ANGELUA); mx[9]=0; mx[10]=cos(dir*ANGELUA); mx[11]=0;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
        break;
    
    case 't':
        mx[0]=1; mx[1]=0; mx[2]=0; mx[3]=0;
        mx[4]=0; mx[5]=1; mx[6]=0; mx[7]=dir*DESPL;
        mx[8]=0; mx[9]=0; mx[10]=1; mx[11]=0;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
        break;

    default:
        break;
    }
}

void form_matrix_z(double *mx, int dir){
    if(analisis_mode){
        mx[0]=1; mx[1]=0; mx[2]=0; mx[3]=0;
        mx[4]=0; mx[5]=1; mx[6]=0; mx[7]=0;
        mx[8]=0; mx[9]=0; mx[10]=1; mx[11]=dir*DESPL;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
    }else{
        switch (aldaketa)
        {
        case 'r':
            mx[0]=cos(dir*ANGELUA); mx[1]=-sin(dir*ANGELUA); mx[2]=0; mx[3]=0;
            mx[4]=sin(dir*ANGELUA); mx[5]=cos(dir*ANGELUA); mx[6]=0; mx[7]=0;
            mx[8]=0; mx[9]=0; mx[10]=1; mx[11]=0;
            mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
            break;
        
        case 't':
            mx[0]=1; mx[1]=0; mx[2]=0; mx[3]=0;
            mx[4]=0; mx[5]=1; mx[6]=0; mx[7]=0;
            mx[8]=0; mx[9]=0; mx[10]=1; mx[11]=dir*DESPL;
            mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
            break;

        default:
            break;
        }
    }
}

void form_resize_matrix(double *mx){
    switch (resize)
    {
    case 0: //Txikitu
        mx[0]=1/SCALE; mx[1]= 0; mx[2]=0; mx[3]= 0;
        mx[4]=0; mx[5]= 1/SCALE; mx[6]= 0; mx[7]= 0;
        mx[8]=0; mx[9]= 0; mx[10]= 1/SCALE; mx[11]= 0;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
        break;
    
    case 1: //Handitu
        mx[0]=SCALE; mx[1]= 0; mx[2]=0; mx[3]= 0;
        mx[4]=0; mx[5]=SCALE; mx[6]= 0; mx[7]= 0;
        mx[8]=0; mx[9]= 0; mx[10]= SCALE; mx[11]= 0;
        mx[12]=0; mx[13]=0; mx[14]=0; mx[15]=1;
        break;

    default:
        break;
    }
}

void form_rodriguez_matrix(double * mx, char rot, int dir){
     double x, y, z, cosa, sina;

    switch (rot){
        case 'y':
            x = camera->m_obj[0];
            y = camera->m_obj[4];
            z = camera->m_obj[8];
            break;

        case 'x':
            x = camera->m_obj[1];
            y = camera->m_obj[5];
            z = camera->m_obj[9];
            break;

        case 'z':
            x = camera->m_obj[2];
            y = camera->m_obj[6];
            z = camera->m_obj[10];
            break;

        default:
            break;
    }

    cosa = cos(dir*ANGELUA);
    sina = sin(dir*ANGELUA);

    mx[0] = (cosa) + ((1-cosa)*(x*x)); mx[1] = ((1-cosa)*(x*y)) - (z*sina); mx[2] = ((1-cosa)*(x*z)) + (y*sina); mx[3]=0;
    mx[4] = ((1-cosa)*(x*y)) + (z*sina); mx[5] = (cosa) + ((1-cosa)*(y*y)); mx[6] = ((1-cosa)*(y*z)) - (x*sina); mx[7]=0;
    mx[8] = ((1-cosa)*(x*z)) - (y*sina); mx[9] = ((1-cosa)*(y*z)) + (x*sina); mx[10] = (cosa) + ((1-cosa)*(z*z)); mx[11] = 0;
    mx[12] = 0; mx[13] = 0; mx[14] = 0; mx[15] = 1;

}

void biraketa_anilisis(int dir, char axis){
    double mt_mat[16], mt_at[16], mr[16], in_between[16];
    punto lau_zut ;
    camera_obj* new_cam = (camera_obj *)malloc(sizeof(camera_obj));

    //Get 4. Zutabea
    lau_zut.x=sel_ptr->mptr->m[3];
    lau_zut.y=sel_ptr->mptr->m[7];
    lau_zut.z=sel_ptr->mptr->m[11];

    //Calculate mt_mat
    mt_mat[0] = 1; mt_mat[1] = 0; mt_mat[2] = 0; mt_mat[3] = -lau_zut.x;
    mt_mat[4] = 0; mt_mat[5] = 1; mt_mat[6] = 0; mt_mat[7] = -lau_zut.y;
    mt_mat[8] = 0; mt_mat[9] = 0; mt_mat[10] = 1; mt_mat[11] = -lau_zut.z;
    mt_mat[12] = 0; mt_mat[13] = 0; mt_mat[14] = 0; mt_mat[15] = 1;

    //Calculate mt_at
    mt_at[0] = 1; mt_at[1] = 0; mt_at[2] = 0; mt_at[3] = lau_zut.x;
    mt_at[4] = 0; mt_at[5] = 1; mt_at[6] = 0; mt_at[7] = lau_zut.y;
    mt_at[8] = 0; mt_at[9] = 0; mt_at[10] = 1; mt_at[11] = lau_zut.z;
    mt_at[12] = 0; mt_at[13] = 0; mt_at[14] = 0; mt_at[15] = 1;

    //Calculate mr
    form_rodriguez_matrix(&(mr[0]),axis,dir);


    //Calculate New Camera
    matrix_calc(&(in_between[0]),&(mr[0]),&(mt_mat[0]));
    matrix_calc(&(in_between[0]),&(mt_at[0]),&(in_between[0]));
    matrix_calc(&(new_cam->m_obj[0]),&(in_between[0]),&(camera->m_obj[0]));

    m_esa_calc(new_cam);

    //Update Camera
    new_cam->prev_camera = camera;
    camera = new_cam;

    printf("\n");
    print_matrizea(camera->m_obj);
}

void x_aldaketa(int dir)
{
    double mx_new[16];
    form_matrix_x(mx_new,dir);
    if (camera_selected==0){
        if(!analisis_mode){
            mlist* result = (mlist *)malloc(sizeof(mlist));
            double * mptr = sel_ptr->mptr->m;//Triangeluaren Matrizearen Lehenengo elemntua
            if(ald_lokala==1){
                matrix_calc(&(result->m[0]),&(mptr[0]),&(mx_new[0]));
            }else{
                matrix_calc(&(result->m[0]),&(mx_new[0]),&(mptr[0]));
            }
            result->hptr= sel_ptr->mptr;
            sel_ptr->mptr= result;
            if (sel_ptr->is_cam && object_as_camera){
                update_camera(&(mx_new[0]));
            }
        }else{
            biraketa_anilisis(dir, 'x');
        }
    }else if (camera_selected && !object_as_camera){
        update_camera(&(mx_new[0]));
    }

}

void y_aldaketa(int dir)
{
    double mx_new[16];
    form_matrix_y(mx_new,dir);
    if (camera_selected==0){
        if (analisis_mode == 0){
            mlist* result = (mlist *)malloc(sizeof(mlist));
            double * mptr = sel_ptr->mptr->m;//Triangeluaren Matrizearen Lehenengo elemntua
            if(ald_lokala==1){
                matrix_calc(&(result->m[0]),&(mptr[0]),&(mx_new[0]));
            }else{
                matrix_calc(&(result->m[0]),&(mx_new[0]),&(mptr[0]));
            }
            result->hptr= sel_ptr->mptr;
            sel_ptr->mptr= result;
            if (sel_ptr->is_cam && object_as_camera){
                update_camera(&(mx_new[0]));
            }
        }else{
            biraketa_anilisis(dir, 'y');
        }
    }else if (camera_selected && !object_as_camera){
        update_camera(&(mx_new[0]));
    }
}

void z_aldaketa(int dir)
{
    double mx_new[16];
    form_matrix_z(mx_new,dir);
    if (camera_selected==0){
        if (!analisis_mode){
            mlist* result = (mlist *)malloc(sizeof(mlist));
            double * mptr = sel_ptr->mptr->m;//Triangeluaren Matrizearen Lehenengo elemntua
            if(ald_lokala==1){
                matrix_calc(&(result->m[0]),&(mptr[0]),&(mx_new[0]));
            }else{
                matrix_calc(&(result->m[0]),&(mx_new[0]),&(mptr[0]));
            }
            result->hptr= sel_ptr->mptr;
            sel_ptr->mptr= result;
            if (sel_ptr->is_cam && object_as_camera){
                update_camera(&(mx_new[0]));
            }
        }else{
            update_camera(&(mx_new[0]));
        }
        printf("\n");
        print_matrizea(camera->m_obj);
    }else if (camera_selected && !object_as_camera){
        update_camera(&(mx_new[0]));
    }
}

void undo()
{  
    if (camera_selected==0)
    {
        if (sel_ptr->mptr!=NULL)sel_ptr->mptr = sel_ptr->mptr->hptr;
    }else{
        if (camera->prev_camera!=NULL) camera = camera->prev_camera;
    }

}

void resize_mx(){
    mlist* result = (mlist *)malloc(sizeof(mlist));
    double mx_new[16];
    double * mptr = sel_ptr->mptr->m;//Triangeluaren Matrizearen Lehenengo elemntua
    form_resize_matrix(mx_new);
    if(ald_lokala==1){
        matrix_calc(&(result->m[0]),&(mptr[0]),&(mx_new[0]));
    }else{
        matrix_calc(&(result->m[0]),&(mx_new[0]),&(mptr[0]));
    }
    result->hptr= sel_ptr->mptr;
    sel_ptr->mptr= result;
}

void look_at_obj(){
    double e_at[3],z_c[3],u_zc[3],x_c[3],y_c[3],v_up[3];
    double e_at_n,u_zc_n,v_up_n;
    punto lau_zut;
    camera_obj* new_cam = (camera_obj *)malloc(sizeof(camera_obj));

    lau_zut.x=sel_ptr->mptr->m[3];
    lau_zut.y=sel_ptr->mptr->m[7];
    lau_zut.z=sel_ptr->mptr->m[11];
    v_up[0]=camera->m_obj[1];
    v_up[1]=camera->m_obj[5];
    v_up[2]=camera->m_obj[9];

    v_up_n = calc_normal(&(v_up[0]));
    normalize(&(v_up[0]),&(v_up[0]),v_up_n);

    //mxp(&lau_zut,modelview,lau_zut);
    //Calculate new camera
    //Calculate Zc
    e_at[0]=camera->m_obj[3]-lau_zut.x;
    e_at[1]=camera->m_obj[7]-lau_zut.y;
    e_at[2]=camera->m_obj[11]-lau_zut.z;

    e_at_n = calc_normal(&(e_at[0]));
    normalize(&(z_c[0]),&(e_at[0]),e_at_n);

    new_cam->m_obj[2]=z_c[0]; 
    new_cam->m_obj[6]=z_c[1]; 
    new_cam->m_obj[10]=z_c[2]; 
    new_cam->m_obj[14]=0;

    //Calculate Xc
    u_zc[0]=(v_up[1]*z_c[2])-(v_up[2]*z_c[1]);
    u_zc[1]= -((v_up[2]*z_c[0])-(v_up[0]*z_c[2])); 
    u_zc[2]= (v_up[0]*z_c[1])-(v_up[1]*z_c[0]);

    u_zc_n = calc_normal(&(u_zc[0]));
    normalize(&(x_c[0]),&(u_zc[0]),u_zc_n);

    new_cam->m_obj[0]=x_c[0]; 
    new_cam->m_obj[4]=x_c[1]; 
    new_cam->m_obj[8]=x_c[2];
    new_cam->m_obj[12]=0;

    //Calculate Yc
    y_c[0]=(z_c[1]*x_c[2])-(z_c[2]*x_c[1]);
    y_c[1]= -((z_c[2]*x_c[0])-(z_c[0]*x_c[2]));
    y_c[2]= (z_c[0]*x_c[1])-(z_c[1]*x_c[0]);
    new_cam->m_obj[1]=y_c[0];
    new_cam->m_obj[5]=y_c[1];
    new_cam->m_obj[9]=y_c[2];
    new_cam->m_obj[13]=0;

    //Set Position
    new_cam->m_obj[3]=camera->m_obj[3];
    new_cam->m_obj[7]=camera->m_obj[7];
    new_cam->m_obj[11]=camera->m_obj[11];
    new_cam->m_obj[15]=1;

    //cam_obj_ptr->mptr->m[3]=camera->m_obj[3];
    //cam_obj_ptr->mptr->m[7]=camera->m_obj[7];
    //cam_obj_ptr->mptr->m[15]=camera->m_obj[15];

    //Calculate ESA
    m_esa_calc(new_cam);

    //Update Camera
    new_cam->prev_camera = camera;
    camera = new_cam;

    print_matrizea(camera->m_obj);
}

// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;
switch(key)
	{
	case 13: 
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_faces) 
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
    case 'c':
        camera_selected = 1-camera_selected;
        break;
    case 'C':
        if (!object_as_camera){
            object_as_camera = 1;
            convert_object_to_camera();
        }else{
            object_as_camera = 0;
            sel_ptr->is_cam = 0;
            camera = camera->prev_camera;
        }
        break;

    case 'G':
        //TODO
        // Hegazkin modua <------> Analisi modua
        analisis_mode = 1-analisis_mode;
        if(analisis_mode){
            look_at_obj();
        }
        update_cam_obj();
        break;

    case 'n':
        show_normal_v = 1-show_normal_v;
        break;

    case 'b':
        back_culling = 1-back_culling;
        //TODO: Atze aurpegiak ezabatu / Marraztu
        break;

    case 'p':
        //TODO: Perspektiba / Paraleloa
        perspective = 1-perspective;
        break;

	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
		break;
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
	case 't':
	        aldaketa = 't';
		break;
	case 'r':
		aldaketa = 'r';
		break;
    case 's':
        resize = 0;
        resize_mx();
        //TODO Txikitu objektua
        break;
    case 'S':
        resize = 1;
        resize_mx();
        //TODO Handitu obketua
        break;
	case 'g':
		if (ald_lokala == 1) ald_lokala = 0;
		    else ald_lokala = 1;
		break;
        case 'x':
                x_aldaketa(1);
                break;
        case 'y':
                y_aldaketa(1);
                break;
        case 'z':
                z_aldaketa(1);
                break;
        case 'X':
                x_aldaketa(-1);
                break;
        case 'Y':
                y_aldaketa(-1);
                break;
        case 'Z':
                z_aldaketa(-1);
                break;
        case 'u':
                undo();
                break;
	case 'f':
	        /*Ask for file*/
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
	        read_from_file(fitxiz);
	        indexx = 0;
                break;
       /* case 'S':  // save to file
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
                if ((obj_file = fopen(fitxiz, "w")) == NULL)
                         {
                         printf("ezin fitxategia ireki\n");
                         }
                     else
                         {
                         for (i =0; i < sel_ptr->num_triangles; i++)
                            {
                            fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                 sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z, 
                                 sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                 sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z, 
                                 sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                 sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z, 
                                 sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                            }
                         fclose(obj_file);
                         }
                break; */
        case 9: /* <TAB> */
            if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                sel_ptr = sel_ptr->hptr;
                /*The selection is circular, thus if we move out of the list we go back to the first element*/
                if (sel_ptr == 0) sel_ptr = foptr->hptr;
                indexx =0; // the selected polygon is the first one
                }
            break;
	case 27:  // <ESC>
		exit( 0 );
		break;
	default:
		printf("%d %c\n", key, key );
	}

// The screen must be drawn to show the new triangle
glutPostRedisplay();
}

void set_camera(){
    double * cam_mx=camera->m_obj;
    cam_mx[0]=1; cam_mx[1]=0; cam_mx[2]=0; cam_mx[3]=0;
    cam_mx[4]=0; cam_mx[5]=1; cam_mx[6]=0; cam_mx[7]=0;
    cam_mx[8]=0; cam_mx[9]=0; cam_mx[10]=1; cam_mx[11]=300;
    cam_mx[12]=0; cam_mx[13]=0; cam_mx[14]=0; cam_mx[15]=1;
    m_esa_calc(camera);
}

void load_projection_mx(){
    float r,t,l,b,n,f;
    r = 5; t = 5; n = 5; l = -5; b = -5, f = 500; 
    projection_mx[0] = (2*n)/(r-l); projection_mx[1] = 0; projection_mx[2] = (r+l)/(r-l); projection_mx[3] = 0;
    projection_mx[4] = 0; projection_mx[5] = (2*n)/(t-b); projection_mx[6] = (t+b)/(t-b); projection_mx[7] = 0;
    projection_mx[8] = 0; projection_mx[9] = 0; projection_mx[10] = -(f+n)/(f-n); projection_mx[11] = (-2*f*n)/(f-n);
    projection_mx[12] = 0; projection_mx[13] = 0; projection_mx[14] = -1; projection_mx[15] = 0;
}

int main(int argc, char** argv)
{
int retval;

	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	glutInitWindowSize ( 500, 500 );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow( "KBG/GO praktika" );
	
	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */ 
        retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
        if (retval !=1) 
            {
            printf("Ez dago texturaren fitxategia (testura.ppm)\n");
            exit(-1);
            }
        
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        denak = 1;
        lineak =1;
        objektuak = 1;
        foptr = 0;
        sel_ptr = 0;
        aldaketa = 'r';
        ald_lokala = 1;
        perspective = 0;
        back_culling = 1;
        camera = (camera_obj *)malloc(sizeof(camera_obj));
        camera_selected=0;
        object_as_camera=0;
        show_normal_v = 0;
        analisis_mode=0;
        set_camera();
        load_projection_mx();
        //matrix_calc(camera_by_mx,camera->m_esa,ca)
        if (argc>1) read_from_file(argv[1]);
            else{
                read_from_file("k.obj");
                printf("Object Loaded");
                //foptr->mptr->m[3] = 250;
                //foptr->mptr->m[11] = 250;
                //read_from_file("z.txt");
                //foptr->mptr->m[3] = -250;
                //read_from_file("z.txt");
                //foptr->mptr->m[7] = 150;
            }
        //read_from_file("cam.txt");
        //cam_obj_ptr = foptr;
        //update_cam_obj();
        //sel_ptr = sel_ptr->hptr;
	glutMainLoop();

	return 0;   
}