/******************************************+
*                                          *
*  CSCI420 - Computer Graphics USC         *
*  Author: Rodrigo Castiel                 *
*                                          *
+*******************************************/
#pragma once

#include "scene_object.h"

/*********************************************************
* This file provides classes for rendering primitive objs
* such as axis, squares, mirrors, spheres/domes, terrain.
* 
* The meshes contain the geometry by itself (vertex info)
* and the inherit from class Mesh.
*
* The corresponding objects store not only the mesh, but
* also material, texture and motion attributes.
* 
* EVERYTHING here was made FROM SCRATCH and took several
* hours!
*
* Rodrigo Castiel, 02/23/2016.
*********************************************************/

// +----------------------------------------------------------+
// | Table of useful objects  | Lighting | Texture | Material |
// | 1. class AxisObject      |    off   |   not   |    no    |
// | 2. class Mirror TODO     |    off   |   yes   |    no    |
// | 3. class GridObject      |    off   |   no    |    no    |
// | 4. class TexturedSphere  |    opt   |   yes   |    no    |
// | 5. class TexturedTerrain |    opt   |   yes   |    no    |
// +----------------------------------------------------------+
// Constructors' default is Object(pipelineProgram, programHandle)

// The basic RGB axis used to show object center and orientation.
class AxisMesh : public Mesh
{
public:
  AxisMesh() { }
  virtual ~AxisMesh() { }

  void Load();
};
//-------------------------------------------------------------
// Square - typically used in mirrors.
class SquareMesh : public Mesh
{
public:
  SquareMesh() { }
  virtual ~SquareMesh() { }

  void Load();
};
//-------------------------------------------------------------
// Rectangular grid defining some plane.
class GridMesh : public Mesh
{
public:
  GridMesh(int w, int h) 
  : mWidth(w), mHeight(h) 
  { }
  
  virtual ~GridMesh() { }
  void Load();

private:
  int mWidth;
  int mHeight;
};
//-------------------------------------------------------------

class TexturedTerrainMesh : public Mesh
{
public:
  TexturedTerrainMesh(int w, int h) 
  : mWidth(w), mHeight(h)
  { }
  
  void Load(ImageIO* heightmap);

  virtual ~TexturedTerrainMesh() { }

  private:
    int mWidth;
    int mHeight;
};
//-------------------------------------------------------------
// Textured sphere - typically used for creating a dome.
class TexturedSphereMesh : public Mesh
{
public:
  TexturedSphereMesh(int detailLevel) 
  : mDetailLevel(detailLevel)
  { }
  
  virtual ~TexturedSphereMesh() { }

  void Load(bool completeDome);

private:
  int mDetailLevel;

};

///////////////////////////////////////////////////////////////////////////////////////////////////

class AxisObject : public SceneObject
{
public:
  AxisObject(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
  : SceneObject(pipelineProgram, programHandle) 
  { }

  void Load()
  {
    AxisMesh* mesh = new AxisMesh();
    mesh->SetProgramHandle(mProgramHandle);
    mesh->Load();
    mMesh = mesh;
    SceneObject::SetScale(5.0f, 5.0f, 5.0f);
    mIsMeshOwner = true;
  }

  virtual ~AxisObject() { }
};

class Mirror : public SceneObject
{
public:
  Mirror(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
  : SceneObject(pipelineProgram, programHandle) 
  { }

  void Load()
  {
    // TexturedSphereMesh* mesh = new TexturedSphereMesh(32);
    // mesh->SetProgramHandle(mProgramHandle);
    // mesh->Load();
    // mesh->UploadGLBuffers();
    // mMesh = mesh;
    // mUsingLighting = false;
    // mTexture = new Texture();
    // mIsMeshOwner = true;
  }

  virtual ~Mirror()
  {
    delete mTexture;
  }
};


class GridObject : public SceneObject
{
public:
  GridObject(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
  : SceneObject(pipelineProgram, programHandle) 
  { }

  void Load(int w = 11, int h = 11)
  {
    GridMesh* mesh = new GridMesh(w, h);
    mesh->SetProgramHandle(mProgramHandle);
    mesh->Load();
    mMesh = mesh;

    SceneObject::SetScale(w, 1.0f, h);
    mIsMeshOwner = true;
  }

  virtual ~GridObject() { }
};

class TexturedSphere : public SceneObject
{
public:
  TexturedSphere(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
  : SceneObject(pipelineProgram, programHandle) 
  { }

  void Load(const std::string& fileName, bool completeDome = true)
  {
    TexturedSphereMesh* mesh = new TexturedSphereMesh(64);
    mesh->SetProgramHandle(mProgramHandle);
    mesh->Load(completeDome);
    mMesh = mesh;

    mTexture = new Texture();
    mTexture->Load(fileName);

    mMaterial = new Material( glm::vec3(0.18), 
                              glm::vec3(0.8), 
                              glm::vec3(0.02) );

    mUsingLighting = true;
    mIsMeshOwner = true;
  }

  virtual ~TexturedSphere()
  {
    delete mTexture;
    delete mMaterial;
  }
};

class TexturedTerrain : public SceneObject
{
public:
  TexturedTerrain(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
  : SceneObject(pipelineProgram, programHandle)
  { }

  void Load(const std::string& heightmapFileName, const std::string& textureFileName)
  {
    ImageIO* source = new ImageIO();
    TexturedTerrainMesh* mesh;
  
    // Try to load the heightmap.
    if (source->loadJPEG(heightmapFileName.c_str()) == ImageIO::OK)
    {
      mesh = new TexturedTerrainMesh(source->getWidth(), source->getHeight());
      mesh->SetProgramHandle(mProgramHandle);
      mesh->Load(source);
    }
    else
    {
      mesh = new TexturedTerrainMesh(20, 20);
      mesh->SetProgramHandle(mProgramHandle);
      mesh->Load(nullptr);
    }

    mTexture = new Texture();
    mTexture->Load(textureFileName);

    mMaterial = new Material( glm::vec3(0.18), 
                              glm::vec3(0.8), 
                              glm::vec3(0.02) );
    mUsingLighting = true;

    mMesh = mesh;
    mIsMeshOwner = true;
  }

  virtual ~TexturedTerrain()
  {
    delete mTexture;
    delete mMaterial;
  }

};
