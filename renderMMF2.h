/**
 * \file renderMMF.h
 * 
 * \brief This file contains the prototype and implementation of strategies
 * to move the foveae.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date April 2019
 *
 * This file is part of projectFoveaCuda software.
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version. This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details. You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RENDERMMF_H
#define RENDERMMF_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef __CUDACC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#ifdef __unix__   // Unix
#include <GL/gl.h>  //GLdouble, Glint, glDrawPixels
#include <GL/glu.h>
#include <GL/glut.h> //gluUnProject
#else // Apple 
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#include <GLUT/glut.h>
#endif
/*#else  // Windows
#include <w32api/GL/gl.h> //GLdouble, Glint, glDrawPixels
#include <w32api/GL/glu.h> 
#include <w32api/GL/glut.h> //gluUnProject
#endif*/


using namespace cv;

/**
 * \struct MMF
 *
 * \brief Struct for fovea using OpenMP and Cuda libraries.
 */
struct RenderMMF{
  
  //
  // Methods
  //
  
  /**
   * \fn Point getDelta( int k, int m, Point w, Point u, Point f )
   *
   * \brief Calculates the initial pixel to build MMF.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) to build the fovea
   *
   * \return Return the initial pixel on the both axis of level k to build MMF.
   */
#ifdef __CUDACC__
  __device__
#endif
  Point getDelta( int k, int m, Point w, Point u, Point f );

  /**
   * \fn Point getSize( int k, int m, Point w, Point u )
   *
   * \brief Calculates the final pixel to build MMF.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *
   * \return Return the final pixel on the both axis of level k to build MMF.
   */
#ifdef __CUDACC__
  __device__
#endif
  Point getSize( int k, int m, Point w, Point u );

  /**
   * \fn std::vector< GLubyte > calcColors( int k, int m, Point d, Point s, Cena* cena, 
   * Luz* luz, Vetor* lookfrom, Vetor* lookat, GLdouble model[16], GLdouble proj[16], 
   * GLint view[4] )
   * 
   * \brief Calculates the colors of the defined region.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        d - Start point of the region to be analyzed
   *        s - End point of the region to be analyzed
   *        cena - Cena que sera aplicado o ray tracing                                                                                                                                        
   *        luz - Luz no objeto                                                                                                                                                                
   *        lookfrom - posicao da camera                                                                                                                                                       
   *        lookat - posicao para onde esta apontada a camera                                                                                                                                  
   *        model, proj, view - matrizes modelview, projection e viewport
   *
   * \return RGB colors to compose the scene
   */
  std::vector< GLubyte > calcColors( int k, int m, Point d, Point s, Cena* cena, Luz* luz, 
				     Vetor* lookfrom, Vetor* lookat, GLdouble model[16], 
				     GLdouble proj[16], GLint view[4] );
  
  /**
   * \fn Mat MMF_CPU( int k, int m, Point w, Point u, Point f,
   * Cena* cena, Luz* luz, Vetor* lookfrom, Vetor* lookat,                                                                                                            
   * GLdouble model[16], GLdouble proj[16], GLint view[4] )
   *
   * \brief Calculates the levels of MMF method using CPUs.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) of the fovea
   *        cena - Cena que sera aplicado o ray tracing                                                                                                                                        
   *        luz - Luz no objeto                                                                                                                                                                
   *        lookfrom - posicao da camera                                                                                                                                                       
   *        lookat - posicao para onde esta apontada a camera                                                                                                                                  
   *        model, proj, view - matrizes modelview, projection e viewport                                                                                                                      
   
   *
   * \return Return the level of MMF method.
   */
  Mat MMF_CPU( int k, int m, Point w, Point u, Point f, Cena* cena, 
	       Luz* luz, Vetor* lookfrom, Vetor* lookat, GLdouble model[16], 
	       GLdouble proj[16], GLint view[4] );
  
  /**
   * \fn Point mapLevel2Image( int k, int m, Point w, Point u, Point f, Point px )
   *
   * \brief Calculates the position of pixel on the level to image.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) of the fovea
   *        px - Pixel (x, y) that we want to map.
   *
   * \return Return the position of pixel on the both axis to image.
   */
  Point mapLevel2Image( int k, int m, Point w, Point u, Point f, Point px );
  
  /**
   * \fn Mat foveated( int m, Point w, Point u, Point f,
   * Cena* cena, Luz* luz, Vetor* lookfrom, Vetor* lookat,                                                                                                            
   * GLdouble model[16], GLdouble proj[16], GLint view[4] )
   *
   * \brief Calculates the levels of MMF method using CPUs.
   *
   * \param m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) of the fovea
   *        cena - Cena que sera aplicado o ray tracing                                                                                                                                        
   *        luz - Luz no objeto                                                                                                                                                                
   *        lookfrom - posicao da camera                                                                                                                                                       
   *        lookat - posicao para onde esta apontada a camera                                                                                                                                  
   *        model, proj, view - matrizes modelview, projection e viewport                                                                                                                      
   *        method - If (0) by default will be considered MMF_CPU, else (1) will be
   *        considered MMF_GPU
   *
   * \return Return the level of MMF method.
   */
  Mat foveated( int m, Point w, Point u, Point f, Cena* cena, 
		Luz* luz, Vetor* lookfrom, Vetor* lookat, GLdouble model[16], 
		GLdouble proj[16], GLint view[4], int method );
  
};

#endif

// Implementation

/**
 * \fn Point getDelta( int k, int m, Point w, Point u, Point f )
 *
 * \brief Calculates the initial pixel to build MMF.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) to build the fovea
 *
 * \return Return the initial pixel on the both axis of level k to build MMF.
 */
#ifdef __CUDACC__
__device__
#endif
Point 
RenderMMF::getDelta( int k, int m, Point w, Point u, Point f ){
  int dx = int( k * ( u.x - w.x + ( 2 * f.x ) ) )/ ( 2 * m );
  int dy = int( k * ( u.y - w.y + ( 2 * f.y ) ) )/ ( 2 * m );
#ifdef DEBUG
  std::cout << "Delta: ( " << dx << ", " << dy << " ) " << std::endl;  
#endif
  return Point( dx, dy );
}

/**
 * \fn Point getSize( int k, int m, Point w, Point u )
 *
 * \brief Calculates the final pixel to build MMF.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *
 * \return Return the final pixel on the both axis of level k to build MMF.
 */
#ifdef __CUDACC__
__device__
#endif
Point 
RenderMMF::getSize( int k, int m, Point w, Point u ){
  int sx = ((m * u.x) + (w.x * k) - (k * u.x)) / m;
  int sy = ((m * u.y) + (w.y * k) - (k * u.y)) / m;
#ifdef DEBUG
  std::cout << "Size: ( " << sx << ", " << sy << " ) " << std::endl;  
#endif
  return Point( sx, sy );
}

/**
 * \fn std::vector< GLubyte > calcColors( int k, int m, Point d, Point s, Cena* cena, 
 * Luz* luz, Vetor* lookfrom, Vetor* lookat, GLdouble model[16], GLdouble proj[16], 
 * GLint view[4] )
 * 
 * \brief Calculates the colors of the defined region.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        d - Start point of the region to be analyzed
 *        s - End point of the region to be analyzed
 *        cena - Cena que sera aplicado o ray tracing                                                                                                                                        
 *        luz - Luz no objeto                                                                                                                                                                
 *        lookfrom - posicao da camera                                                                                                                                                       
 *        lookat - posicao para onde esta apontada a camera                                                                                                                                  
 *        model, proj, view - matrizes modelview, projection e viewport
 *
 * \return RGB colors to compose the scene
 */
std::vector< GLubyte > 
RenderMMF::calcColors( int k, int m, Point d, Point s, Cena* cena, Luz* luz, 
		       Vetor* lookfrom, Vetor* lookat, GLdouble model[16], 
		       GLdouble proj[16], GLint view[4] ){
  std::vector< GLubyte > colors (3);
  colors[0] = 0.0; colors[1] = 0.0; colors[2] = 0.0;
  Cena* cena_auxiliar = new Cena();
  Objeto* objeto_salvo = new Objeto();
  double t_aux;
  double t_auxMax = -1.0;
  // Loop parameters
  int divisionsRegion = 2*( m - (k+1) ) + 1;
  Point gapBetweenSteps = Point( (int)((s.x + d.x)/ divisionsRegion), (int)((s.y + d.y)/ divisionsRegion) );
  Point stepByStepRegion = Point( (int)(gapBetweenSteps.x / 2), (int)(gapBetweenSteps.y / 2) ); 
  //std::cout << "d: " << "(" << d.x << ", " << d.y << ")" << std::endl;
  //std::cout << "s: " << "(" << d.x + s.x << ", " << d.y + s.y << ")" << std::endl;
  //std::cout << "divisionsRegion: " << divisionsRegion << std::endl;
  //std::cout << "GapBetweenSteps: " << "(" << gapBetweenSteps.x << ", " << gapBetweenSteps.y << ")" << std::endl;
  //std::cout << "stepByStepRegion: " << "(" << stepByStepRegion.x << ", " << stepByStepRegion.y << ")" << std::endl;
  
  // Scene
  //for ( int sx = d.x; sx <= d.x + 1; sx++ ){ // Getting the first pixel of region
  //  for ( int sy = d.y; sy < d.y + 1; sy++ ){
      
  for ( int sx = d.x; sx <= d.x + s.x; sx++ ){ // Getting all pixels of region
    for ( int sy = d.y; sy < d.y + s.y; sy++ ){
  
  // Inteligent sampling, but wasn't finished
  //for ( int sx = d.x + stepByStepRegion.x; sx < s.x + d.x; sx+=gapBetweenSteps.x ){
  //  for ( int sy = d.y + stepByStepRegion.y; sy < s.y + d.y; sy+=gapBetweenSteps.y ){
      t_aux = -1.0; // Absurd value t
      //Looking for lookat's
      GLdouble x, y, z;
      GLint realy = view[3] - (GLint)sy;
      //Read the window z value from the z-buffer
      glReadPixels( sx, realy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z );
      gluUnProject((GLdouble)sx, (GLdouble)realy, 1, model, proj, view, &x, &y, &z);
      Vetor* lookatAux = new Vetor();
      lookatAux->valores_vetor((double)x, (double)y, (double)z);
      
      Raio* r = new Raio(); // Ray
      
      int tamanho = cena->size_objetos_pilha();
      for (int k = 0; k < tamanho; k++){ // Loop scene
	Objeto* obj = cena->excluir_objetos_pilha();
	r->atualizar_vetores(lookfrom, lookatAux, obj);
	double t = r -> calcula_t();
	if (t > 0.0){         //Verificar se o t encontrado e maior que zero                                                                                                                   
	  if (t_aux == -1.0){ //Nao tem nenhum valor maior que 0 que seja menor que -1                                                                                                         
	    t_aux = t;
	    objeto_salvo = obj;
	  }
	  else{
	    if (t < t_aux){           //Verificar se e menor que os t's encontrados                                                                                                            
	      t_aux = t;
	      objeto_salvo = obj;
	    }
	  }
	}
	cena_auxiliar->incluir_objetos_pilha(obj);
      }
      
      Cena* tmp = cena;
      cena = cena_auxiliar;
      cena_auxiliar = tmp;

      //if ( t_aux != -1.0 )
      //std::cout << t_aux << std::endl;

      if (t_aux > 0.0){ // Sphere has been intercepted
	r->atualizar_vetores(lookfrom, lookatAux, objeto_salvo);
	Vetor* int_esfera = new Vetor();
	int_esfera = r->interseccao_esfera(t_aux);
	luz->atualizar_vetores_auxiliares(int_esfera, objeto_salvo->posicao_esfera(), lookfrom);
	Vetor* cores_objeto = new Vetor();
	// Fixed texture
	cores_objeto = objeto_salvo->cor_esfera();
	// Random texture
	//Limites da cor do objeto
	//Vetor* range_areia_superior = new Vetor();  //Vetor de limite superior da cor                                                                                                        
	//range_areia_superior->valores_vetor(223.0, 246.0, 143.0);                                                                                                                            
	//Vetor* range_areia_inferior = new Vetor();  //Vetor de limite inferior da cor                                                                                                        
	//range_areia_inferior->valores_vetor(139.0, 129.0, 76.0);                                                                                                                             
	//Textura* cor_areia = new Textura(range_areia_superior, range_areia_inferior);                                                                                                        
	//Aplicando a textura ao objeto                                                                                                                                                        
	//objeto_salvo->modificar_cor_pixel(cor_areia);                                                                                                                                        
	//cores_objeto = objeto_salvo->cor_esfera();
	
	Vetor* cor_luz = new Vetor();
	double valor_luz_vermelha = luz->calcula_luz_red();
	double valor_luz_verde = luz->calcula_luz_green();
	double valor_luz_azul = luz->calcula_luz_blue();
	
	cor_luz->valores_vetor(valor_luz_vermelha, valor_luz_verde, valor_luz_azul);
	
	if ( t_auxMax < t_aux ){
	  t_auxMax = t_aux;
	  colors[2] = (GLubyte)(cor_luz->vx() * (cores_objeto->vx()/cores_objeto->norma())) ;
	  colors[1] = (GLubyte)(cor_luz->vy() * (cores_objeto->vy()/cores_objeto->norma())) ;
	  colors[0] = (GLubyte)(cor_luz->vz() * (cores_objeto->vz()/cores_objeto->norma())) ;
	}

      }
      
    }
  }
  
  return colors;
}


/**
 * \fn Mat MMF_CPU( int k, int m, Point w, Point u, Point f,
 * Cena* cena, Luz* luz, Vetor* lookfrom, Vetor* lookat,                                                                                                            
 * GLdouble model[16], GLdouble proj[16], GLint view[4] )
 *
 * \brief Calculates the levels of MMF method using CPUs.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) of the fovea
 *        cena - Cena que sera aplicado o ray tracing                                                                                                                                        
 *        luz - Luz no objeto                                                                                                                                                                
 *        lookfrom - posicao da camera                                                                                                                                                       
 *        lookat - posicao para onde esta apontada a camera                                                                                                                                  
 *        model, proj, view - matrizes modelview, projection e viewport                                                                                                                      
 *
 *
 * \return Return the level of MMF method.
 */
Mat
RenderMMF::MMF_CPU( int k, int m, Point w, Point u, Point f, Cena* cena, 
		    Luz* luz, Vetor* lookfrom, Vetor* lookat, GLdouble model[16], 
		    GLdouble proj[16], GLint view[4] ){
  // Checking conditions
  assert( ( w.x > 0 ) && ( w.x < u.x ) );
  assert( ( w.y > 0 ) && ( w.y < u.y ) );
  assert( ( u.x > 0 ) && ( u.y > 0 ) );
  assert( m >= 1 );
  
  // Creating img
  //GLubyte image[w.x][w.y][3];
  Mat image(w.x, w.y, CV_8UC3);
  
  // Delta e Size
  Point d = getDelta( k, m, w, u, f );
  Point s = getSize( k, m, w, u );
  
  Point regions = Point( (int)(s.x/w.x), (int)(s.y/w.y) ); 
  //Point regions = Point( (int)(s.x - d.x)/w.x, (int)(s.y - d.y)/w.y );
#ifdef DEBUG
  std::cout << "regions: " << "(" << regions.x << ", " << regions.y << ")" << std::endl;
#endif
  
  std::vector< GLubyte > colors;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/ 
#endif
  for ( int wi = 0; wi < w.x; wi++ ){ // Completing W
    for ( int wj = 0; wj < w.y; wj++ ){
      Point startRegion = Point(d.x + wi*regions.x, d.y + wj*regions.y );
      Point finishRegion = Point( regions.x, regions.y );
      //std::cout << "startRegion: " << "(" << startRegion.x << ", " << startRegion.y << ")" << std::endl;
      //std::cout << "finishRegion: " << "(" << finishRegion.x << ", " << finishRegion.y << ")" << std::endl;
      colors = calcColors( k, m, startRegion, finishRegion, cena, luz, lookfrom, lookat, model, proj, view );
      Vec3b rgb = Vec3b((uchar)colors[0], (uchar)colors[1], (uchar)colors[2]);
      image.at<Vec3b>((w.x-1) - wi, wj) = rgb;
      colors.clear();      
    }
  }

  return image;

}

/**
 * \fn Point mapLevel2Image( int k, int m, Point w, Point u, Point f, Point px )
 *
 * \brief Calculates the position of pixel on the level to image.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) of the fovea
 *        px - Pixel (x, y) that we want to map.
 *
 * \return Return the position of pixel on the both axis to image.
 */
Point 
RenderMMF::mapLevel2Image( int k, int m, Point w, Point u, Point f, Point px ){
  int _px = ( (k * w.x) * (u.x - w.x) + (2 * k * w.x * f.x) + (2 * px.x) * ( (m * u.x) - (k * u.x) + (k * w.x) ) )/ (2 * m * w.x);
  int _py = ( (k * w.y) * (u.y - w.y) + (2 * k * w.y * f.y) + (2 * px.y) * ( (m * u.y) - (k * u.y) + (k * w.y) ) )/ (2 * m * w.y);
#ifdef DEBUG
  std::cout << "Map: ( " << _px << ", " << _py << " ) " << std::endl;  
#endif
  return Point( _px, _py );
}

/**
 * \fn Mat foveated( int m, Point w, Point u, Point f,
 * Cena* cena, Luz* luz, Vetor* lookfrom, Vetor* lookat,                                                                                                            
 * GLdouble model[16], GLdouble proj[16], GLint view[4] )
 *
 * \brief Calculates the levels of MMF method using CPUs.
 *
 * \param m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) of the fovea
 *        cena - Cena que sera aplicado o ray tracing                                                                                                                                        
 *        luz - Luz no objeto                                                                                                                                                                
 *        lookfrom - posicao da camera                                                                                                                                                       
 *        lookat - posicao para onde esta apontada a camera                                                                                                                                  
 *        model, proj, view - matrizes modelview, projection e viewport                                                                                                                      
 *        method - If (0) by default will be considered MMF_CPU, else (1) will be
 *        considered MMF_GPU
 *
 * \return Return the level of MMF method.
 */
Mat 
RenderMMF::foveated( int m, Point w, Point u, Point f, Cena* cena, 
		     Luz* luz, Vetor* lookfrom, Vetor* lookat, GLdouble model[16], 
		     GLdouble proj[16], GLint view[4], int method ){
  Mat imgFoveated(u.x, u.y, CV_8UC3);
  if ( method == 0 ){ // MMF_CPU
#ifdef _OPENMP
#pragma omp parallel for schedule(static, m+1) // Schedule(static, m+1) keeps the order
#endif
    for ( int k = 0; k < m + 1; k++ ){ // Levels
      Mat imgLevel = MMF_CPU( k, m, w, u, f, cena, luz, lookfrom, lookat, model, proj, view );
      // Mapping levels to foveated image
      Point initial = mapLevel2Image( k, m, w, u, f, Point( 0, 0 ) );
      Point final = mapLevel2Image( k, m, w, u, f, Point( w.x, w.y ) );
#ifdef DEBUG
      std::cout << "(xi, yi) = (" << initial.x << ", " << initial.y << ")" << std::endl;
      std::cout << "(xf, yf) = (" << final.x << ", " << final.y << ")" << std::endl;
#endif
      Rect roi = Rect( initial.x, initial.y, final.x - initial.x, final.y - initial.y );
      if ( k < m ){ // Copying levels to foveated image
        resize( imgLevel, imgLevel, Size(final.x - initial.x, final.y - initial.y), 0, 0, CV_INTER_LINEAR );
        imgLevel.copyTo( imgFoveated( roi ) );
      }
      else
        imgLevel.copyTo( imgFoveated( roi ) );
      // Paint levels
      cv::rectangle(imgFoveated, cv::Point(initial.x, initial.y), cv::Point(final.x - 1, final.y - 1), cv::Scalar(255, 255, 255));
    }
  }

  return imgFoveated;
}
