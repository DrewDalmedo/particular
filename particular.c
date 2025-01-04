#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <stdlib.h>
#include <math.h>

#define SCREEN_WIDTH    736
#define SCREEN_HEIGHT   584
#define MAX_PARTICLES   1000000   
#define G               0.1       // Gravitational constant (adjusted for simulation scale)
#define PULL_STRENGTH   5.0       // Strength of the pull force
#define ANCHOR_MASS     1000.0    // Mass of the anchor point
#define ANCHOR_RADIUS   1.5       // Increased radius for initial particle positions
#define TIME_STEP       0.075     // Reduced time step for better accuracy

#define BG_COLOR        WHITE
#define FG_COLOR        BLACK

typedef struct {
  Uint8 r;
  Uint8 g;
  Uint8 b;
} Color;

const Color WHITE = {
  255,
  255,
  255
};
    
const Color BLACK = {
  0,
  0,
  0
};

typedef struct {
    double x, y;        // Position
    double vx, vy;      // Velocity
    double ax, ay;      // Acceleration
    double mass;

    int active;         // whether particle is active

    int anchored;       // whether particle is anchor
    double anchor_x, anchor_y; // Anchor position
    double anchor_mass; // Mass of the anchor point

    Color color;        // Particle color
} Particle;

// Convert an SDL_Surface to grayscale
SDL_Surface* ConvertToGrayscale(SDL_Surface* surface) {

    // Convert source surface to RGBA32 format
    SDL_Surface* formatted_surface = SDL_ConvertSurfaceFormat(surface, SDL_PIXELFORMAT_RGBA32, 0);
    if (!formatted_surface) {
        printf("Failed to convert surface format: %s\n", SDL_GetError());
        return NULL;
    }

    SDL_Surface* grayscale = SDL_CreateRGBSurfaceWithFormat(0, formatted_surface->w, formatted_surface->h, 8, SDL_PIXELFORMAT_INDEX8);
    if (!grayscale) {
        printf("Failed to create grayscale surface: %s\n", SDL_GetError());
        SDL_FreeSurface(formatted_surface);
        return NULL;
    }

    // Set grayscale palette
    if (grayscale->format->palette) {
        SDL_Color palette[256];
        for (int i = 0; i < 256; i++) {
            palette[i].r = palette[i].g = palette[i].b = i;
            palette[i].a = 255;
        }
        SDL_SetPaletteColors(grayscale->format->palette, palette, 0, 256);
    } else {
        printf("Failed to get grayscale palette.\n");
        SDL_FreeSurface(grayscale);
        SDL_FreeSurface(formatted_surface);
        return NULL;
    }

    // Copy pixels to grayscale surface
    Uint32* src_pixels = (Uint32*)formatted_surface->pixels;
    Uint8* dst_pixels = (Uint8*)grayscale->pixels;
    int src_pitch = formatted_surface->pitch / 4; // Pixels per row
    int dst_pitch = grayscale->pitch; // Bytes per row (since each pixel is 1 byte)

    for (int y = 0; y < formatted_surface->h; y++) {
        for (int x = 0; x < formatted_surface->w; x++) {
            Uint32 pixel = src_pixels[y * src_pitch + x];
            Uint8 r, g, b, a;
            SDL_GetRGBA(pixel, formatted_surface->format, &r, &g, &b, &a);

            // Calculate luminance
            Uint8 gray = (Uint8)(0.3 * r + 0.59 * g + 0.11 * b);
            dst_pixels[y * dst_pitch + x] = gray;
        }
    }

    SDL_FreeSurface(formatted_surface);
    return grayscale;
}

// Simple edge detection a la Sobel operator
Uint8* EdgeDetection(SDL_Surface* grayscale_surface) {
    int width = grayscale_surface->w;
    int height = grayscale_surface->h;
    Uint8* edge_data = (Uint8*)malloc(width * height);
    if (!edge_data) {
        printf("Failed to allocate memory for edge data.\n");
        return NULL;
    }

    Uint8* pixels = (Uint8*)grayscale_surface->pixels;
    int pitch = grayscale_surface->pitch; // Bytes per row

    // Initialize edge_data to zero
    memset(edge_data, 0, width * height);

    // Sobel operator kernels
    int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    int Gy[3][3] = {
        {-1, -2, -1},
        { 0,  0,  0},
        { 1,  2,  1}
    };

    // Edge detection
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int sumX = 0;
            int sumY = 0;

            // Apply kernels
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    int pixel = pixels[(y + ky) * pitch + (x + kx)];
                    sumX += pixel * Gx[ky + 1][kx + 1];
                    sumY += pixel * Gy[ky + 1][kx + 1];
                }
            }

            // Calculate gradient magnitude
            int magnitude = abs(sumX) + abs(sumY);
            if (magnitude > 255) magnitude = 255;

            edge_data[y * width + x] = (Uint8)magnitude;
        }
    }

    return edge_data;
}

int main(int argc, char* argv[]) {
    Color bg_color = BG_COLOR;
    Color fg_color = FG_COLOR;

    // Check for image path argument
    if (argc < 2) {
        printf("Usage: %s <image_path>\n", argv[0]);
        return 1;
    }
    char* image_path = argv[1];

    // Initialize SDL2 and SDL2_image
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("SDL Initialization Failed: %s\n", SDL_GetError());
        return 1;
    }

    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG)) {
        printf("SDL_image Initialization Failed: %s\n", IMG_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("Particle Image Renderer",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        SCREEN_WIDTH, SCREEN_HEIGHT, 0);
    if (!window) {
        printf("Window Creation Failed: %s\n", SDL_GetError());
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
    if (!renderer) {
        printf("Renderer Creation Failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    // Enable blending mode for transparency
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    // Load the image
    SDL_Surface* image_surface = IMG_Load(image_path);
    if (!image_surface) {
        printf("Image Loading Failed: %s\n", IMG_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    // Convert image to grayscale
    SDL_Surface* grayscale_surface = ConvertToGrayscale(image_surface);
    SDL_FreeSurface(image_surface); // Free the original surface
    if (!grayscale_surface) {
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    // Perform edge detection
    Uint8* edge_data = EdgeDetection(grayscale_surface);
    if (!edge_data) {
        SDL_FreeSurface(grayscale_surface);
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    // Get image dimensions
    int img_width = grayscale_surface->w;
    int img_height = grayscale_surface->h;

    // Allocate particles array on the heap
    Particle* particles = (Particle*)malloc(sizeof(Particle) * MAX_PARTICLES);
    // If we couldn't allocate memory for particles, exit safely
    if (!particles) {
        printf("Failed to allocate memory for particles.\n");
        free(edge_data);
        SDL_FreeSurface(grayscale_surface);
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    // Initialize particles to zero
    memset(particles, 0, sizeof(Particle) * MAX_PARTICLES);

    int particle_count = 0;

    // Scale factors to fit the image to the screen
    double scale_x = (double)SCREEN_WIDTH / img_width;
    double scale_y = (double)SCREEN_HEIGHT / img_height;

    int max_particles_reached = 0;

    // Threshold for edge detection (adjust as needed)
    Uint8 edge_threshold = 80;

    // Create particles at edge positions
    for (int y = 0; y < img_height && !max_particles_reached; y++) {
        for (int x = 0; x < img_width; x++) {
            if (particle_count >= MAX_PARTICLES) {
                max_particles_reached = 1;
                break;
            }

            Uint8 magnitude = edge_data[y * img_width + x];

            if (magnitude >= edge_threshold) {
                double pos_x = x * scale_x;
                double pos_y = y * scale_y;

                // Initialize particle
                particles[particle_count].mass = 1.0;
                particles[particle_count].active = 1;
                particles[particle_count].anchored = 1;
                particles[particle_count].anchor_x = pos_x;
                particles[particle_count].anchor_y = pos_y;
                particles[particle_count].anchor_mass = ANCHOR_MASS;

                // Set radius for initial distance from anchor
                double radius = ANCHOR_RADIUS;
                double angle = ((double)(rand() % 360)) * M_PI / 180.0;

                // Position particle at a distance 'radius' from anchor point
                particles[particle_count].x = pos_x + radius * cos(angle);
                particles[particle_count].y = pos_y + radius * sin(angle);

                // Calculate speed for circular orbit
                double speed = sqrt(G * particles[particle_count].anchor_mass / radius);

                // Set initial velocity perpendicular to the radius vector
                particles[particle_count].vx = -speed * sin(angle);
                particles[particle_count].vy = speed * cos(angle);

                particles[particle_count].ax = 0;
                particles[particle_count].ay = 0;

                particle_count++;
            }
        }
    }

    printf("Loaded %d particles from image edges.\n", particle_count);

    // Free resources
    free(edge_data);
    SDL_FreeSurface(grayscale_surface);

    int running = 1;
    SDL_Event event;
    int mouse_x, mouse_y;
    Uint32 mouse_state;
    int action = 1;  // Default action

    while (running) {
        // Event handling
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT)
                running = 0;
            else if (event.type == SDL_KEYDOWN) {
                // Select action with number keys 1-9
                if (event.key.keysym.sym >= SDLK_1 && event.key.keysym.sym <= SDLK_9) {
                    action = event.key.keysym.sym - SDLK_0;
                    printf("Action changed to %d\n", action);
                }
            }
        }

        // Get mouse state
        mouse_state = SDL_GetMouseState(&mouse_x, &mouse_y);

        // Left-click to perform action
        if (mouse_state & SDL_BUTTON(SDL_BUTTON_LEFT)) {
            if (particle_count < MAX_PARTICLES) {
                if (action == 1) {
                    // Action 1: Spawn particles at cursor
                    particles[particle_count].x = mouse_x;
                    particles[particle_count].y = mouse_y;
                    particles[particle_count].vx = ((double)(rand() % 200) - 100) / 50.0;
                    particles[particle_count].vy = ((double)(rand() % 200) - 100) / 50.0;
                    particles[particle_count].mass = 1.0;
                    particles[particle_count].active = 1;
                    particles[particle_count].anchored = 0;
                    particles[particle_count].anchor_mass = 0.0;
                    particle_count++;
                } else if (action == 2) {
                    // Action 2: Spawn anchored particle that orbits around the cursor position
                    particles[particle_count].anchor_x = mouse_x;
                    particles[particle_count].anchor_y = mouse_y;
                    particles[particle_count].anchor_mass = ANCHOR_MASS;

                    // Set initial position offset from anchor point
                    double angle = ((double)(rand() % 360)) * M_PI / 180.0;
                    double radius = ANCHOR_RADIUS;  // Distance from anchor point (can be varied)
                    particles[particle_count].x = mouse_x + radius * cos(angle);
                    particles[particle_count].y = mouse_y + radius * sin(angle);

                    // Calculate velocity for circular orbit
                    double speed = sqrt(G * particles[particle_count].anchor_mass / radius);
                    particles[particle_count].vx = -speed * sin(angle);
                    particles[particle_count].vy = speed * cos(angle);

                    particles[particle_count].mass = 1.0;
                    particles[particle_count].active = 1;
                    particles[particle_count].anchored = 1;
                    particle_count++;
                }
                // Additional actions can be added here
            }
        }

        // Clear screen
        
        SDL_SetRenderDrawColor(renderer, bg_color.r, bg_color.g, bg_color.b, 255);
        SDL_RenderClear(renderer);

        // Reset accelerations
        for (int i = 0; i < particle_count; i++) {
            if (!particles[i].active) continue;
            particles[i].ax = 0;
            particles[i].ay = 0;
        }

        // Calculate gravitational forces
        for (int i = 0; i < particle_count; i++) {
            if (!particles[i].active) continue;

            // For anchored particles, add gravitational pull towards anchor point
            if (particles[i].anchored) {
                double dx = particles[i].anchor_x - particles[i].x;
                double dy = particles[i].anchor_y - particles[i].y;
                double dist_sq = dx * dx + dy * dy;
                if (dist_sq < 1e-6) dist_sq = 1e-6;  // Prevent division by zero or extremely large forces
                double dist = sqrt(dist_sq);
                double force = G * particles[i].anchor_mass * particles[i].mass / dist_sq;
                particles[i].ax += force * dx / dist / particles[i].mass;
                particles[i].ay += force * dy / dist / particles[i].mass;
            }
            else {
                // Gravitational interaction with other particles
                for (int j = i + 1; j < particle_count; j++) {
                    if (!particles[j].active) continue;
                    if (particles[j].anchored) continue; // Skip anchored particles
                    double dx = particles[j].x - particles[i].x;
                    double dy = particles[j].y - particles[i].y;
                    double dist_sq = dx * dx + dy * dy + 0.01;  // Avoid division by zero
                    double dist = sqrt(dist_sq);
                    double force = G * particles[i].mass * particles[j].mass / dist_sq;
                    double fx = force * dx / dist;
                    double fy = force * dy / dist;

                    particles[i].ax += fx / particles[i].mass;
                    particles[i].ay += fy / particles[i].mass;
                    particles[j].ax -= fx / particles[j].mass;
                    particles[j].ay -= fy / particles[j].mass;
                }
            }
        }

        // Apply pulling force towards cursor when right-click is held
        if (mouse_state & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
            for (int i = 0; i < particle_count; i++) {
                if (!particles[i].active) continue;
                double dx = mouse_x - particles[i].x;
                double dy = mouse_y - particles[i].y;
                double dist_sq = dx * dx + dy * dy + 0.01;
                double dist = sqrt(dist_sq);
                double force = PULL_STRENGTH * particles[i].mass / dist_sq;
                particles[i].ax += force * dx / dist;
                particles[i].ay += force * dy / dist;
            }
        }

        // Update positions and velocities using fixed time step
        for (int i = 0; i < particle_count; i++) {
            if (!particles[i].active) continue;
            particles[i].vx += particles[i].ax * TIME_STEP;
            particles[i].vy += particles[i].ay * TIME_STEP;
            particles[i].x += particles[i].vx * TIME_STEP;
            particles[i].y += particles[i].vy * TIME_STEP;

            // Keep particles within the screen bounds
            if (particles[i].x < 0) particles[i].x = 0;
            if (particles[i].x >= SCREEN_WIDTH) particles[i].x = SCREEN_WIDTH - 1;
            if (particles[i].y < 0) particles[i].y = 0;
            if (particles[i].y >= SCREEN_HEIGHT) particles[i].y = SCREEN_HEIGHT - 1;
        }

        // Draw particles
        SDL_SetRenderDrawColor(renderer, fg_color.r, fg_color.g, fg_color.b, 255);
        for (int i = 0; i < particle_count; i++) {
            if (!particles[i].active) continue;
            SDL_RenderDrawPoint(renderer, (int)particles[i].x, (int)particles[i].y);
        }

        SDL_RenderPresent(renderer);
        SDL_Delay(16);  // ~60 FPS
    }

    // Clean up and exit
    free(particles);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    IMG_Quit();
    SDL_Quit();
    return 0;
}
