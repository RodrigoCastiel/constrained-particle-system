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
  
  GridObject* originGrid = new GridObject(pipelineProgram, programHandle);

  TexturedSphere* sky = new TexturedSphere(pipelineProgram, programHandle);
  sky->SetScale(100.0, 100.0, 100.0);

  Camera* mainCamera = new Camera(pipelineProgram, programHandle);
  // mainCamera->Scale(2, 2, );
  mCameras.push_back(mainCamera);

  sky->SetRotVelocity(2e-6, 1e-6, -1e-5f);

  originGrid->Load(15, 15);
  originGrid->SetPosition(0, -1, 0);

  sky->Load("textures/outer_space.jpg");
  sky->SetLighting(false);

  mObjects.push_back(originGrid);
  mObjects.push_back(sky);

  // Initialie particle system here.
  mParticleSystem = new ParticleSystem(pipelineProgram, programHandle);
  mParticleSystem->Setup(11, {0, -.5}, {.0, .0});

  mInitialized = true;
}

void CPSScene::Render()
{
  Scene::Render();
  mParticleSystem->Render();
}

void CPSScene::Animate()
{
  Scene::Animate();
  mParticleSystem->Animate();
}

void CPSScene::OnMouseLeftClick(int x, int y, int w, int h, int state)
{
  Camera* camera = mCameras[mCurrentCamera];
  glm::vec3 r = camera->ComputeRayAt(x, y, w, h);
  glm::vec3 C = camera->GetCenterCoordinates();
  mParticleSystem->SelectParticle(C, r, state);
}

void CPSScene::OnMouseLeftDrag(int x, int y, int w, int h)
{
  Camera* camera = mCameras[mCurrentCamera];
  glm::vec3 r = camera->ComputeRayAt(x, y, w, h);
  glm::vec3 C = camera->GetCenterCoordinates();
  mParticleSystem->DragParticle(C, r);
}

void CPSScene::OnMouseRightClick(int x, int y, int w, int h)
{

}

void CPSScene::Clean()
{
  Scene::Clean();
  delete mParticleSystem;
}
