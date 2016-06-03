#include "scene.h"

void Scene::Init(BasicPipelineProgram* pipelineProgram, GLuint programHandle)
{
  mPipelineProgram = pipelineProgram;
  mProgramHandle = programHandle;
  mInitialized = true;
}

void Scene::Render()
{
  if (Scene::IsInitialized())
  { 
    // Position current camera.
    mCameras[mCurrentCamera]->Position();
    OpenGLMatrix& currentView = mCameras[mCurrentCamera]->GetViewMatrix();

    // Set all lights.
    for (auto light : mLights)
    {
      light->Position(currentView);
    }

    // Render all objects.
    for (auto object : mObjects)
    {
      object->Render();
    }
  }
}

void Scene::ChangeCamera()
{
  int n = mCameras.size();
  mCurrentCamera = (mCurrentCamera + 1) % n;
}

void Scene::Animate()
{
  for (auto camera : mCameras)
  {
    camera->Animate();
  }

  for (auto light : mLights)
  {
    light->Animate();
  }

  for (auto object : mObjects)
  {
    object->Animate();
  }
}

void Scene::ReshapeScreen(int w, int h)
{
  glViewport(0, 0, w, h);

  // Update each camera's projection.
  for (auto camera : mCameras)
  {
    camera->Project(0, 0, w, h);
  }
}

void Scene::Clean()
{
  for (auto object : mObjects)
  {
    delete object;
  }

  for (auto light : mLights)
  {
    delete light;
  }

  for (auto camera : mCameras)
  {
    delete camera;
  }
}