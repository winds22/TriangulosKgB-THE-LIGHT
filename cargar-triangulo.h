
typedef struct punto
{
    double x, y, z, u,v;
} punto;

typedef struct normal_v
{
    float x, y, z;
}normal_v;

typedef struct hiruki
{
    punto p1,p2,p3;
    normal_v bek_normal;
} hiruki;

int cargar_triangulos(char *fitxiz, int *hkopptr, hiruki **hptrptr);


int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr, unsigned char **rgbptr);