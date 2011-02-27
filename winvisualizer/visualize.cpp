#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <gl/gl.h>
#include <gl/glu.h>
#include <math.h>
#include <vector>

#pragma warning(disable:4996)

#define DEFAULT_FILENAME "sample.txt"
#define CLASSNAME "Particle Simulation"
#define eps 0.1
#define SCALE 200
#define FPS 15
#define MIN_SIZE 100

struct particle_t { float x, y; };

//
//  timer
//
double read_timer( )
{	
    static bool initialized = false;
    static double dfreq;
    static LARGE_INTEGER seconds0;
    LARGE_INTEGER temp;
    if( !initialized )
    {
        QueryPerformanceFrequency (&temp);
        dfreq = 1.0/temp.QuadPart;
        QueryPerformanceCounter (&seconds0);
        initialized = true;
    }
    QueryPerformanceCounter (&temp);
    return (temp.QuadPart - seconds0.QuadPart)*dfreq;
}	

//
//  window message handler
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
    if( msg == WM_CLOSE || ( msg == WM_KEYDOWN && wParam == VK_ESCAPE ) )
    {
        PostQuitMessage (0);
        return 0;
    }
	
    if( msg == WM_SIZE )
    {
        int window_sx = LOWORD (lParam);
        int window_sy = HIWORD (lParam);
        glViewport (0, 0, window_sx, window_sy);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        return 0;
    }

    return DefWindowProc( hWnd, msg, wParam, lParam );
}

int WINAPI WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int iCmdShow )
{
    char *filename = lpCmdLine[0] ? lpCmdLine : DEFAULT_FILENAME;

    FILE *f = fopen( filename, "r" );
    if( f == NULL )
        return 1;
    
    int n;
    float size;
    fscanf( f, "%d%g", &n, &size );

    particle_t p;
    std::vector<particle_t> particles;
    while( fscanf( f, "%g%g", &p.x, &p.y ) == 2 )
        particles.push_back( p );
    fclose( f );

    int nframes = particles.size( ) / n;
    if( nframes == 0 )
        return 2;
    
    int window_size = (int)((size+2*eps)*SCALE);
    window_size = window_size > MIN_SIZE ? window_size : MIN_SIZE;
    
    //
    //  create window
    //
    WNDCLASS wc;
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = WndProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wc.hCursor = LoadCursor (NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH) GetStockObject (BLACK_BRUSH);
    wc.lpszMenuName = NULL;
    wc.lpszClassName = CLASSNAME;
    RegisterClass( &wc );

    HWND hWindow = CreateWindow( CLASSNAME, CLASSNAME,
        WS_CAPTION | WS_POPUPWINDOW | WS_VISIBLE | WS_THICKFRAME | WS_SYSMENU | WS_MAXIMIZEBOX,
        0, 0, window_size+8, window_size+34, NULL, NULL, hInstance, NULL );

    HDC hDC = GetDC (hWindow);

    //
    //  init OpenGL
    //
    PIXELFORMATDESCRIPTOR pfd;
    ZeroMemory( &pfd, sizeof (pfd) );
    pfd.nSize = sizeof( pfd );
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    pfd.cDepthBits = 16;
    pfd.cAlphaBits = 16;
    pfd.iLayerType = PFD_MAIN_PLANE;
    int format = ChoosePixelFormat( hDC, &pfd );
    SetPixelFormat( hDC, format, &pfd );
    HGLRC hRC = wglCreateContext( hDC );
    wglMakeCurrent( hDC, hRC );

    glEnable( GL_POINT_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

    glPointSize( 1 );
    glClearColor( 1, 1, 1, 1 );

    //
    //  runtime loop
    //
    MSG msg;
    while( 1 )
    {
        if( PeekMessage( &msg, NULL, 0, 0, PM_REMOVE ) )
        {
            if (msg.message == WM_QUIT) 
                break;
            else 
            {
                TranslateMessage( &msg );
                DispatchMessage( &msg );
            }
        } 
        else 
        {
            glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
			
            glMatrixMode( GL_PROJECTION );
            glLoadIdentity( );
            gluOrtho2D( -eps, size+eps, -eps, size+eps );
            
            glColor3f( 0.75, 0.75, 0.75 );
            glBegin( GL_LINE_LOOP );
            glVertex2d( 0, 0 );
            glVertex2d( size, 0 );
            glVertex2d( size, size );
            glVertex2d( 0, size );
            glEnd( );
            
            int iframe = (int)(read_timer()*FPS) % nframes;
            particle_t *p = &particles[iframe*n];
        
            glColor3f( 0, 0, 0 );
            glBegin( GL_POINTS );
            for( int i = 0; i < n; i++ )
                glVertex2fv( &p[i].x );
            glEnd( );
			
            glFinish( );
            SwapBuffers( hDC );
        }
    }

    //shutdown
    wglMakeCurrent( NULL, NULL );
    wglDeleteContext( hRC );
    ReleaseDC( hWindow, hDC );
    DestroyWindow( hWindow );
    
    return msg.wParam;
}
