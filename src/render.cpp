#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "texture.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

         if (features.enableTextureMapping && hitInfo.material.kdTexture && !features.enableShading) {
            if (features.extra.enableTransparency && rayDepth < 5) {
                 Ray r;

                 r.origin = ray.origin + 0.000025f * r.direction;
                 r.direction = ray.direction;


                 bvh.intersect(r, hitInfo, features);

                 glm::vec3 firstColor = Lo;
                 glm::vec3 secondColor = getFinalColor(scene, bvh, r, features, rayDepth + 1);

                 Lo = hitInfo.material.transparency * firstColor + secondColor * (1 - hitInfo.material.transparency);
                 drawRay(r, Lo);
            }
            return acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
        }

        if (enableDebugDraw && features.extra.enableMotionBlur) {
            Ray r = ray;
            glm::vec3 orgn = r.origin;
            HitInfo current = hitInfo;
            for (int i = -20; i <= 20; i++) {
                r.t = std::numeric_limits<float>::max();
                bvh.intersect(r, hitInfo, features);
                r.origin = orgn + (float)i * features.vector / 10000.0f;
                drawRay(r, { 0, 0, 1 });
            }
            hitInfo = current;
        }

        if (features.enableRecursive) {
            if (rayDepth < 3) {
                if (hitInfo.material.ks[0] > 0.0f || hitInfo.material.ks[1] > 0.0f || hitInfo.material.ks[2] > 0.0f) {
                    Lo += hitInfo.material.ks * getFinalColor(scene, bvh, computeReflectionRay(ray, hitInfo), features, rayDepth + 1); 
                }
            }
        }

        drawRay(ray, Lo);



        if (rayDepth < 5 && features.extra.enableTransparency) {
            Ray r = Ray { ray.origin + ray.direction * ray.t, ray.direction };
            r.origin = r.origin + 0.000025f * r.direction;
            glm::vec3 result = Lo * hitInfo.material.transparency;
            return result + (1 - hitInfo.material.transparency) * getFinalColor(scene, bvh, r, features, rayDepth + 1);
        }

        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

glm::vec3 boxFilter(Screen& source, std::vector<glm::vec3>& pixels, int i, int j, int filterSize)
{
    glm::ivec2 windowResolution = source.resolution();
    glm::vec3 sum = glm::vec3(0.0f);
    for (int x = -filterSize; x < filterSize + 1; ++x) {
        if (i + x >= 0 && i + x < windowResolution.x) {
            for (int y = -filterSize; y < filterSize + 1; ++y) {
                if (j + y >= 0 && j + y < windowResolution.y) {
                    sum += pixels[source.indexAt(i + x, j + y)];
                }
            }
        }
    }
        
    sum /= (2 * filterSize + 1) * (2 * filterSize + 1);
    return sum;
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            if (features.extra.enableMotionBlur) {
                glm::vec3 z = glm::vec3(0);
                for (int i = 0; i < 50; i++) {
                    int random = rand() % 41 - 20;
                    glm::vec3 adjustment = features.vector * (float)random / 100.0f;
                    Ray r = camera.generateRay(normalizedPixelPos);
                    r.origin = r.origin + adjustment;

                    z = z + getFinalColor(scene, bvh, r, features);
                }
                z /= 50;
                screen.setPixel(x, y, z);
                continue;
            }
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }
    if (features.extra.enableBloomEffect) {
        float threshold = features.extra.bloomThreshold;
        int filterSize = features.extra.boxSize;
        Screen thresholded = Screen(screen.resolution());
        std::vector<glm::vec3> pixels = screen.pixels();
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x < windowResolution.x; x++) {
                int ind = screen.indexAt(x, y);
                glm::vec3 pixel = pixels[ind];
                if (pixel[0] < threshold) {
                    pixel[0] = 0.0f;
                }
                if (pixel[1] < threshold) {
                    pixel[1] = 0.0f;
                }
                if (pixel[2] < threshold) {
                    pixel[2] = 0.0f;
                }
                thresholded.setPixel(x, y, pixel);
            }
        }
        std::vector<glm::vec3> thresholdedPixels = thresholded.pixels();
        if (features.extra.onlyBloomDisplay) {
            for (int y = 0; y < windowResolution.y; y++) {
                for (int x = 0; x < windowResolution.x; x++) {
                    screen.setPixel(x, y, boxFilter(thresholded, thresholdedPixels, x, y, filterSize));
                }
            }
        } else {
            for (int y = 0; y < windowResolution.y; y++) {
                for (int x = 0; x < windowResolution.x; x++) {
                    int ind = screen.indexAt(x, y);
                    glm::vec3 pixel = pixels[ind];
                    screen.setPixel(x, y, pixel + boxFilter(thresholded, thresholdedPixels, x, y, filterSize));
                }
            }
        }
    }
}