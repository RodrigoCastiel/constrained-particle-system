#include "cps_scene.h"

#define GLM_SWIZZLE
#include <glm/glm.hpp>

void CPSScene::Init(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
{
  Scene::Init(pipelineProgram, programHandle);

  // Set OpenGL Environment.
  glPointSize(2.0);
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  // Default light.
  mLights.push_back(new Light(pipelineProgram, programHandle));
  
  AxisObject* originAxis = new AxisObject(pipelineProgram, programHandle);
  GridObject* originGrid = new GridObject(pipelineProgram, programHandle);

  // Extra.
  TexturedSphere* sky     = new TexturedSphere(pipelineProgram, programHandle);
  TexturedSphere* earth   = new TexturedSphere(pipelineProgram, programHandle);
  TexturedTerrain* terrain = new TexturedTerrain(pipelineProgram, programHandle);

  sky->SetScale(1000.0, 1000.0, 1000.0);
  terrain->SetScale(2200.0f, 1.0f, 2200.0f);
  terrain->SetPosition(0, -20, 0);

  earth->SetPosition(0.0f, 0.0f, 0.0f);
  earth->SetRotation(M_PI/6, M_PI, 0);

  mCameras.push_back(new Camera(pipelineProgram, programHandle));
  mCameras.push_back(new Camera(pipelineProgram, programHandle));

  sky->SetRotVelocity(2e-5, 1e-5, -1e-4f);

  originAxis->Load();
  originGrid->Load(15, 15);
  earth->Load("textures/earth.jpg");
  sky->Load("textures/outer_space.jpg");

  terrain->SetLighting(true);
  sky->SetLighting(false);
  earth->SetLighting(false);

  mObjects.push_back(originAxis);
  mObjects.push_back(originGrid);
  mObjects.push_back(earth);
  mObjects.push_back(sky);

  mInitialized = true;
}

void CPSScene::Render()
{
  Scene::Render();

  // // Render cameras as axis.
  // for (auto camera : mCameras)
  // {
  //   OpenGLMatrix& M = mObjects[0]->GetModelMatrix();
  //   OpenGLMatrix& V = camera->GetViewMatrix();
  //   M = V;
  //   M.Invert();
  //   //M.Scale(10, 10, 10);
  //   mObjects[0]->Render();
  // }
}

void CPSScene::Animate()
{
  Scene::Animate();
}

void CPSScene::OnMouseLeftClick(int x, int y, int w, int h)
{
  Camera* camera = mCameras[mCurrentCamera];
  glm::vec3 r = camera->ComputeRayAt(x, y, w, h);
  glm::vec3 C = camera->GetCenterCoordinates();

  // TODO: Interaction with chain.
}

void CPSScene::OnMouseRightClick(int x, int y, int w, int h)
{

}

void CPSScene::Clean()
{
  Scene::Clean();


}

