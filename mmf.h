/**
 * \file mmf.h
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

#ifndef MMF_H
#define MMF_H

#include <stdio.h>
#include <iostream>
//#include <cuda.h>
//#include <cuda_runtime.h>

using namespace cv;

/**
 * \struct MMFCuda
 *
 * \brief Struct for fovea using cuda library.
 */
struct MMFCuda{
  
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
  Point getSize( int k, int m, Point w, Point u );

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
   * \fn Mat MMF_CPU( Mat img, int k, int m, Point w, Point u, Point f )
   *
   * \brief Calculates the levels of MMF method using CPUs.
   *
   * \param img - Image to be foveated.
   *        k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) of the fovea
   *
   * \return Return the level of MMF method.
   */
  Mat MMF_CPU( Mat img, int k, int m, Point w, Point u, Point f );
  
  /**
   * \fn Mat mmf( Mat img, int k, int m, Point w, Point u, Point f )
   *
   * \brief Calculates the levels of MMF method.
   *
   * \param img - Image to be foveated.
   *        k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) of the fovea
   *
   * \return Return the level of MMF method.
   */
  Mat mmf( Mat img, int k, int m, Point w, Point u, Point f );
  
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
Point 
MMFCuda::getDelta( int k, int m, Point w, Point u, Point f ){
  int dx = int( k * ( u.x - w.x + ( 2 * f.x ) ) )/ ( 2 * m );
  int dy = int( k * ( u.y - w.y + ( 2 * f.y ) ) )/ ( 2 * m );
#ifndef DEBUG
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
Point 
MMFCuda::getSize( int k, int m, Point w, Point u ){
  int sx = ((m * u.x) + (w.x * k) - (k * u.x)) / m;
  int sy = ((m * u.y) + (w.y * k) - (k * u.y)) / m;
#ifndef DEBUG
  std::cout << "Size: ( " << sx << ", " << sy << " ) " << std::endl;  
#endif
  return Point( sx, sy );
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
MMFCuda::mapLevel2Image( int k, int m, Point w, Point u, Point f, Point px ){
  int _px = ( (k * w.x) * (u.x - w.x) + (2 * k * w.x * f.x) + (2 * px.x) * ( (m * u.x) - (k * u.x) + (k * w.x) ) )/ (2 * m * w.x);
  int _py = ( (k * w.y) * (u.y - w.y) + (2 * k * w.y * f.y) + (2 * px.y) * ( (m * u.y) - (k * u.y) + (k * w.y) ) )/ (2 * m * w.y);
#ifndef DEBUG
  std::cout << "Map: ( " << _px << ", " << _py << " ) " << std::endl;  
#endif
  return Point( _px, _py );
}

/**
 * \fn Mat MMF_CPU( Mat img, int k, int m, Point w, Point u, Point f )
 *
 * \brief Calculates the levels of MMF method using CPUs.
 *
 * \param img - Image to be foveated.
 *        k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) of the fovea
 *
 * \return Return the level of MMF method.
 */
Mat 
MMFCuda::MMF_CPU( Mat img, int k, int m, Point w, Point u, Point f ){
  Point d = getDelta( k,  m, w, u, f );
  Point s = getSize( k, m, w, u );
  Mat imgLevel = img( Rect( d.x, d.y, s.x, s.y ) ); // Getting ROI of image
  Mat imgLevelResult;
  if ( k < m )
    resize( imgLevel, imgLevelResult, Size(w.x, w.y), 0, 0, CV_INTER_LINEAR ); // Read page 171 of book Hands on GPU Accelerated Computer Vision With OpenCV And Cuda
#ifdef DEBUG
  imshow( "levels", imgLevelResult );
  waitKey( 0 ); // Waiting enter
#endif
  return imgLevel;
}
  
  
/**
 * \fn Mat mmf( Mat img, int k, int m, Point w, Point u, Point f )
 *
 * \brief Calculates the levels of MMF method.
 *
 * \param img - Image to be foveated.
 *        k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) of the fovea
 *
 * \return Return the level of MMF method.
 */
Mat mmf( Mat img, int k, int m, Point w, Point u, Point f ){
  /*Point d = getDelta( k,  m, w, u, f );
  Point s = getSize( k, m, w, u );
  //Mat imgLevel( s.y - d.y, s.x - d.x, img.type() );
  Mat h_imgLevel = img( Rect( d.x, d.y, s.x, s.y ) ); // Getting ROI of image
  cv::cuda::GpuMat d_imgLevel, d_imgLevelResult; // Declaring images on device
  d_imgLevel.upload( h_imgLevel ); // Uploading imgLevel to device
  if ( k < m )
  cv::cuda::resize( d_imgLevel, d_imgLevelResult, w, cv::INTER_LINEAR); // Read page 171 of book Hands on GPU Accelerated Computer Vision With OpenCV And Cuda*/
  Mat h_imgLevelResult;
  /*d_imgLevelResult.download( h_imgLevelResult );
#ifdef DEBUG
  imshow( "levels", h_imgLevelResult );
  waitKey( 0 ); // Waiting enter
  #endif*/
  return h_imgLevelResult;
}
