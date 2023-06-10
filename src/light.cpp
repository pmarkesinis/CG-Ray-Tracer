#include "light.h"
#include "config.h"
# include <iostream>
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <random>


// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color, float place)
{

    glm::vec3 segDir = segmentLight.endpoint0 - segmentLight.endpoint1;
    position = place * segDir + segmentLight.endpoint1;
    color = place * segmentLight.color0 + (1 - place) * segmentLight.color1;
}


// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color, float placeX, float placeY)
{
    glm::vec3 rowDir = parallelogramLight.edge01;
    glm::vec3 colDir = parallelogramLight.edge02;

    position = (parallelogramLight.v0 + rowDir * placeX) + (colDir * placeY);
    color = (parallelogramLight.color0 * (1 - placeX) + parallelogramLight.color1 * placeX) * (1 - placeY) + (parallelogramLight.color2 * (1 - placeX) + parallelogramLight.color3 * placeX) * (placeY);
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise

/**
* For point light sources you should compute hard shadows. 
* At each intersection point cast a shadow ray towards each light source to determine if the point is in shadow or not. 
* If the point is in shadow for a light source, then that light source should not contribute to the shading computation.
*/
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)

{
        float visibility;

        Ray shadowRay = Ray();
       
        // calculate the offset
        glm::vec3 offset = hitInfo.normal * 0.00001f;
        shadowRay.origin = ray.t * ray.direction + ray.origin + offset;
        shadowRay.direction = glm::normalize(samplePos - shadowRay.origin);
        shadowRay.t = glm::length(samplePos - shadowRay.origin);
        
        // checking if point is in shadow for a light source 
        if (bvh.intersect(shadowRay, hitInfo, features) == true) {
            // in shadow 
            visibility = 0.0f;
        }
        else {
            // light contributes to shading
            visibility = 1.0f;
        }

        // VISUAL DEBUG

        if (visibility == 0.0f) {
            // hidden blue
          drawRay(shadowRay, glm::vec3(0.0, 0.0, 1.0));
     }

      if (visibility == 1.0f) {
            // visible red 
         drawRay(shadowRay, debugColor);
          
      }
        return visibility;
    }

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 color{ 0.0f, 0.0f, 0.0f };
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);

                if (features.enableHardShadow == true) {
                    float visibility = testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo);
                    // returns black if not visible, color if visible 
                    color += visibility * computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                  
                }
                else {
                    color += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                }

              }

            else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);

                if (features.enableSoftShadow == true) {

                    glm::vec3 newPosition; 
                    glm::vec3 newColor;
                    float visibility;
                    int count = 0;
                    float place = 0.0f;

            

                    // RANDOMISER 
                    std::random_device random;  
                    std::default_random_engine genr(random());
                    std::uniform_real_distribution<> distr(0.0, 1.0);

                    // find the direction of the segment
                    glm::vec3 segDir = segmentLight.endpoint0 - segmentLight.endpoint1;
                    // length of the line segment 
                    float segSize = glm::length(segDir) * 10.0f;
                
      
                    for (int x = 0; x < segSize - 1; x++) {
                        count++;
                        place = (float(x) + distr(genr))/ segSize; 
                  
                        // returns position and color for the point lights
                        sampleSegmentLight(segmentLight, newPosition, newColor, place);
                        
                       if (features.enableHardShadow) {
                           visibility = testVisibilityLightSample(newPosition, newColor, bvh, features, ray, hitInfo);
                       } else {
                            visibility = 1.0f;
                       }
                       color += visibility * computeShading(newPosition, newColor, features, ray, hitInfo);
                    } 
                     color = color / float(count);
                }


            }
            else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                if (features.enableSoftShadow == true) {
                    glm::vec3 newPosition;
                    glm::vec3 newColor;
                    float visibility;
                    int count = 0;

                    glm::vec3 rowDir = parallelogramLight.edge01;
                    glm::vec3 colDir = parallelogramLight.edge02;
                     
                    float rowSize = glm::length(parallelogramLight.edge01) * 20.0f;
                    float colSize = glm::length(parallelogramLight.edge02) * 20.0f;


                    // RANDOMISER 
                    std::random_device random;
                    std::default_random_engine genr(random());
                    std::uniform_real_distribution<> distr(0.0, 1.0);
                   
                    for (int x = 0; x < rowSize; x++) {
                        float placeX = (float(x) + distr(genr)) / rowSize;
                        for (int y = 0; y < colSize; y++) {
                            float placeY = (float(y) + distr(genr)) / colSize;
                            count++;
                            sampleParallelogramLight(parallelogramLight, newPosition, newColor, placeX, placeY);

                            if (features.enableHardShadow) {
                                visibility = testVisibilityLightSample(newPosition, newColor, bvh, features, ray, hitInfo);
                            }
                            else {
                                visibility = 1.0f;
                            }
                            color += visibility * computeShading(newPosition, newColor, features, ray, hitInfo);
                        }
                    }
                    color = color / (float) count;
                }
            }
        }
    
        // TODO: replace this by your own implementation of shading
         return color;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
