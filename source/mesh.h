/******************************************+
*                                          *
*  CSCI420 - Computer Graphics USC         *
*  Author: Rodrigo Castiel                 *
*                                          *
+*******************************************/

#pragma once

#include <cmath>

#include <vector>
#include <string>
#include <iostream>

#include "openGLHeader.h"
#include "glutHeader.h"

#include "openGLMatrix.h"
#include "basicPipelineProgram.h"


//  +-------------------------------------------------+
//  |  This class implements a generic 3d mesh to     |
//  |  be rendered using OpenGL core profile.         |
//  |  Its vertices are tightly packed:               |
//  |  [(x y z) (r g b) (nx ny nz) (tx ty)]           |
//  |  This is not related to materials or texture    |
//  |                                                 |
//  |  It is optimized by storing only what it uses.  |
//  |  That is, if it doesn't have texture, the tex-  |
//  |  ture coordinates aren't allocated. The same    |
//  |  with colors and normals.                       |
//  |  It also uses element arrays and the draw mode  |
//  |  is configurable.                               |
//  +-------------------------------------------------+
 
class Mesh
{
public:
  // Constructors.
  Mesh() { }
  //Mesh(const ObjMesh* objMesh);

  // TODO: COMMENT!!
  virtual void Render() const;                    // Renders the geometry at the current origin.
  virtual void Update();                          // Sends the geometry again to the graphics card.

  virtual void UploadGLBuffers();                 // Sends the geometry to the graphics card.

  // Preallocates the geometry arrays as described on the top. It is meant to be used with generic shapes.
  // You  may want to use PositionAt, ColorAt, NormalAt, TexCoordAt and IndexAt to edit the geometry.
  virtual void Preallocate(int numVertices, int numIndices, bool hasColors, bool hasNormals, bool hasTexCoord);

  virtual void LoadFromObj(const std::string& objFilename);

  // Following methods help to load primitive objects such as a curve, a sequence of points with
  // uniform color and probably without lighting.
  virtual void LoadFromPositions(const std::vector<glm::vec3>& positions, GLenum drawMode = GL_POINTS, 
                                 float r = 1, float g = 1, float b = 1);
  virtual void UpdateFromPositions(const std::vector<glm::vec3>& positions);

  // Query methods.
  inline bool IsInitialized() const { return mInitialized; }
  inline bool HasTexCoord()   const { return mHasTexCoord; }
  inline bool HasNormals()    const { return mHasNormals;  }
  inline bool HasColors()     const { return mHasColors;   }

  // Getter and setters.
  inline int GetNumVertices() const { return mNumVertices; }
  inline int GetNumIndices()  const { return mNumIndices;  }
  inline int GetVertexSize()  const { return mVertexSize;  } 

  inline void SetGLDrawMode(GLenum mode) { mDrawMode = mode; };
  inline void SetProgramHandle(GLuint programHandle) { mProgramHandle = programHandle; }

  // Geometry access: don't attempt to access if the corresponding query methods return false.
  GLfloat* PositionAt(int index);   // Address to vertex[index] 3D position.
  GLfloat* ColorAt(int index);      // Address to vertex[index] RGB color.
  GLfloat* NormalAt(int index);     // Address to vertex[index] 3D normal vector.
  GLfloat* TexCoordAt(int index);   // Address to vertex[index] 2D texture coordinates.
  GLuint*  IndexAt(int index);      // Address to element_array[index].

  // Destructor.
  virtual ~Mesh();

protected:
  // Provides a link to the shaders.
  GLuint mProgramHandle;

  // True if it has any geometry loaded in it.
  bool mInitialized { false };

  // Tells geometry information.
  bool mHasTexCoord { false };
  bool mHasNormals  { false };
  bool mHasColors   { false };

  // The vertices, normals, texture and colors array and element array.
  GLfloat* mVertices { nullptr };
  GLuint*  mIndices  { nullptr };
  int mNumVertices { 0 };
  int mNumIndices  { 0 };
  int mVertexSize  { 0 };

  // OpenGL Buffers parameters.
  GLenum mDrawMode { GL_TRIANGLES };
  GLuint mEab { 0 };
  GLuint mVao { 0 };  
  GLuint mVbo { 0 };
}; // Mesh.

inline 
GLfloat* Mesh::PositionAt(int index)
{
  return &mVertices[mVertexSize * index];
}

inline 
GLfloat* Mesh::ColorAt(int index)
{
  return &mVertices[mVertexSize * index + 3];
}

inline 
GLfloat* Mesh::NormalAt(int index)
{
  return &mVertices[mVertexSize*index + (3 + 3*mHasColors)];
}

inline 
GLfloat* Mesh::TexCoordAt(int index)
{
  return &mVertices[mVertexSize*index + (3 + 3*mHasColors + 3*mHasNormals)];
}

inline
GLuint* Mesh::IndexAt(int index)
{
  return &mIndices[index];
}
