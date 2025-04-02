#!/usr/bin/env python3
import os
import imageio.v2 as imageio  # use imageio.v2 for compatibility

def create_video_from_images(image_folder, output_video="output_video.mp4", fps=2):
    # Get a sorted list of image file paths (assumes names sort in the desired order)
    image_files = sorted([
        os.path.join(image_folder, f)
        for f in os.listdir(image_folder)
        if f.lower().endswith((".png", ".jpg", ".jpeg"))
    ])
    
    if not image_files:
        print("No image files found in", image_folder)
        return
    
    print("Found the following images:")
    for img in image_files:
        print(img)
    
    # Load images into frames list
    frames = [imageio.imread(fname) for fname in image_files]
    
    # Save the frames as a video using ffmpeg format
    imageio.mimsave(output_video, frames, fps=fps, format="ffmpeg")
    print("Video saved as", output_video)

if __name__ == "__main__":
    # Change this folder to where your images are stored
    image_folder = "./cc3d_output"  # e.g., "./images"
    create_video_from_images(image_folder, output_video="rotating_graph.mp4", fps=1)
