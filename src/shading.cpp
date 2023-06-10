#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    float normalLightDot = glm::dot(hitInfo.normal, glm::normalize(lightPosition - (ray.origin + ray.t * ray.direction)));
    float normalViewDot = glm::dot(hitInfo.normal, -ray.direction * ray.t);
    if (normalLightDot < 0.0f) {
        return glm::vec3(0, 0, 0);
    }
    glm::vec3 diffuse = lightColor * hitInfo.material.kd * normalLightDot;
    if (normalViewDot < 0.0f) {
        return diffuse;
    }
    
    // update Kd value if texture mapping occurs 
    if (features.enableTextureMapping == true) {
        if (hitInfo.material.kdTexture) {
            Image texelImage = *hitInfo.material.kdTexture.get();
            glm::vec3 updatedKd = acquireTexel(texelImage, hitInfo.texCoord, features);
            // update for diffuse 
            hitInfo.material.kd = updatedKd;
            diffuse = lightColor * hitInfo.material.kd * normalLightDot;
        }
    }

    glm::vec3 v = glm::normalize(-ray.direction * ray.t);
    glm::vec3 r = glm::normalize(glm::reflect(ray.origin + ray.t * ray.direction - lightPosition, hitInfo.normal));
    glm::vec3 specular = lightColor * hitInfo.material.ks * glm::pow(std::abs(glm::dot(v, r)), hitInfo.material.shininess);
    return diffuse + specular;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    glm::vec3 normal = hitInfo.normal;
    if (glm::dot(ray.direction, hitInfo.normal) < 0.0f) {
        normal = -normal;
    }
    glm::vec3 newDirection = ray.direction - 2 * glm::dot(ray.direction, normal) * normal;
    glm::vec3 newOrigin = ray.origin + ray.t * ray.direction - 10 * std::numeric_limits<float>::epsilon() * normal;
    Ray reflectionRay { newOrigin, newDirection, std::numeric_limits<float>::max() };
    return reflectionRay;
}