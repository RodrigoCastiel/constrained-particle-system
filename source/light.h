/*******************************************
*
*  CSCI420 - Computer Graphics USC
*  Author: Rodrigo Castiel
*
*  This class implements the Phong Lighing
*  model.
*
********************************************/
#pragma once

#include "openGLHeader.h"
#include "glutHeader.h"

#include "openGLMatrix.h"
#include "basicPipelineProgram.h"

class Light
{
public:
  Light(BasicPipelineProgram* pipelineProgram, GLuint programHandle) 
  : mPipelineProgram(pipelineProgram), mProgramHandle(programHandle)
  { } 

  // Sets the light in the environment in camera coordinates (sends to shader uniform).
  void Position(OpenGLMatrix& viewMatrix);

  // Updates (for an animation).
  void Animate();

  // SetPipelineProgramParam must be called before rendering.
  inline void SetPipelineProgramParam(BasicPipelineProgram *pipelineProgram, GLuint programHandle) 
  { 
    mProgramHandle   = programHandle;
    mPipelineProgram = pipelineProgram;
  }

private:
  glm::vec3 mPos {0.0f,  200.0f, 0.0f};  // Center in world coordinates.
  glm::vec3 mDir {0.0f, 1.0f,  0.0f};    // Direction vector in world coordinates.

  GLfloat mLa {0.1f};  // Ambient component.
  GLfloat mLd {0.8f};  // Diffuse component.
  GLfloat mLs {0.1f};  // Specular component.

  BasicPipelineProgram *mPipelineProgram {nullptr};
  GLuint mProgramHandle {0};
};