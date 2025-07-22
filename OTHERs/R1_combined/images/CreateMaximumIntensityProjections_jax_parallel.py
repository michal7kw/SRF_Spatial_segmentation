import tifffile
import numpy as np
from pathlib import Path
import concurrent.futures
import jax
import jax.numpy as jnp
import os
import time

# Optional: Configure JAX. By default, JAX will try to use a GPU if available.
# To force CPU: jax.config.update("jax_platform_name", "cpu")
# To ensure float64 precision if needed (usually not for images):
# jax.config.update("jax_enable_x64", True)

def print_jax_devices():
    """Prints available JAX devices."""
    try:
        print("Available JAX devices:")
        for device in jax.devices():
            print(f"  - {device} (Platform: {device.platform}, Device ID: {device.id})")
        # Check for GPU specifically
        gpus = [d for d in jax.devices() if d.platform.lower() == 'gpu' or 'gpu' in str(d.device_kind).lower()] # Check common platform names
        if gpus:
            print(f"JAX GPU(s) detected: {gpus}")
        else:
            print("JAX GPU not detected. JAX will use CPU or other available accelerators.")
    except Exception as e:
        print(f"Could not query JAX devices: {e}")

def process_channel(channel_name, base_path_str, num_z_slices=7):
    """
    Processes a single channel: reads its Z-stack, computes max projection using JAX,
    and saves the result.
    """
    base_path = Path(base_path_str)
    process_start_time = time.time()
    print(f"Starting processing for channel: {channel_name}")
    
    z_stack_images = []
    for z in range(num_z_slices):
        file_path = base_path / f'mosaic_{channel_name}_z{z}.tif'
        try:
            img = tifffile.imread(file_path)
            z_stack_images.append(img)
        except FileNotFoundError:
            print(f"Warning: File not found {file_path}, skipping for channel {channel_name}")
            # If a slice is missing, we might not want to proceed with this channel
            return None 
        except Exception as e:
            print(f"Error reading {file_path} for channel {channel_name}: {e}")
            return None

    if not z_stack_images:
        print(f"No images loaded for channel {channel_name}, skipping max projection.")
        return None

    max_projection_np = None
    try:
        # Convert list of numpy arrays to a JAX array
        # jnp.stack creates a new dimension (axis 0) from the list of 2D arrays
        print(f"  Channel {channel_name}: Stacking {len(z_stack_images)} images for JAX processing...")
        jax_z_stack = jnp.stack(z_stack_images, axis=0)
        
        print(f"  Channel {channel_name}: Computing max projection with JAX...")
        # Compute max projection using JAX
        max_projection_jax = jnp.max(jax_z_stack, axis=0)
        
        # Ensure computation is done and convert back to NumPy array for saving
        max_projection_jax.block_until_ready() # Ensures JAX computation is complete
        max_projection_np = np.asarray(max_projection_jax)
        print(f"  Channel {channel_name}: JAX max projection successful.")

    except Exception as e:
        print(f"Error during JAX processing for channel {channel_name}: {e}")
        print(f"  Channel {channel_name}: Falling back to NumPy for max projection.")
        try:
            if z_stack_images: # Ensure list is not empty
                np_z_stack = np.stack(z_stack_images, axis=0)
                max_projection_np = np.max(np_z_stack, axis=0)
                print(f"  Channel {channel_name}: NumPy max projection successful.")
            else:
                return None # Should not happen if we passed the earlier check
        except Exception as np_e:
            print(f"Error during NumPy fallback for channel {channel_name}: {np_e}")
            return None
            
    if max_projection_np is None:
        print(f"  Channel {channel_name}: Max projection could not be computed.")
        return None

    output_filename = base_path / f'mosaic_{channel_name}_maxproj_jax.tif'
    try:
        print(f"  Channel {channel_name}: Saving max projection to {output_filename}...")
        tifffile.imwrite(output_filename, max_projection_np)
        duration = time.time() - process_start_time
        print(f"Successfully processed and saved: {output_filename} (took {duration:.2f}s)")
        return str(output_filename)
    except Exception as e:
        print(f"Error writing {output_filename} for channel {channel_name}: {e}")
        return None

def main():
    script_start_time = time.time()
    
    print_jax_devices() # Print available JAX devices at the start

    base_path = Path('.')  # Assumes script is run from the images directory
    
    # List of channels to process
    channels = ['DAPI', 'Pcp4', 'PolyT', 'Prox1', 'Snap25', 'TAG', 'Ttr', 'hSETBP1']
    num_z_slices = 7 # As per original script

    # Determine number of workers for parallel processing
    # os.cpu_count() can be None, so handle that.
    cpu_cores = os.cpu_count()
    num_workers = max(1, cpu_cores - 1 if cpu_cores and cpu_cores > 1 else 1)
    # You can cap the number of workers if desired, e.g., num_workers = min(num_workers, 4)
    print(f"Using {num_workers} worker processes.")

    results = []
    # Using ProcessPoolExecutor for CPU-bound tasks (JAX on CPU) and I/O-bound tasks
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks: process_channel function for each channel
        # Pass base_path as a string because Path objects might not always be
        # reliably pickled across processes in all Python/OS environments.
        future_to_channel = {
            executor.submit(process_channel, channel, str(base_path), num_z_slices): channel 
            for channel in channels
        }
        
        for future in concurrent.futures.as_completed(future_to_channel):
            channel = future_to_channel[future]
            try:
                result_path = future.result()
                if result_path:
                    results.append(result_path)
                    # Detailed print is now inside process_channel
                else:
                    print(f"Channel {channel} processing failed or produced no output.")
            except Exception as exc:
                print(f'Channel {channel} generated an exception during future.result(): {exc}')
    
    print("\n--- Processing Summary ---")
    if results:
        print("Successfully generated max projections for:")
        for r_path in results:
            print(f"  - {r_path}")
    else:
        print("No max projections were successfully generated.")
    
    total_duration = time.time() - script_start_time
    print(f"\nTotal script execution time: {total_duration:.2f} seconds.")

if __name__ == '__main__':
    # Set the multiprocessing start method to 'spawn' to avoid issues with JAX
    # and multithreading when using os.fork() (the default on Unix).
    # This should be done before any multiprocessing-related objects are created.
    import multiprocessing
    try:
        multiprocessing.set_start_method('spawn', force=True)
    except RuntimeError:
        # Might have already been set, or context already started.
        # If 'spawn' is already the method, this is fine.
        # If another method was set and context started, this might not take effect.
        # Best practice is to set it once at the very start.
        pass # Or log a warning if strict control is needed

    # This check is crucial for multiprocessing to work correctly on
    # platforms like Windows, preventing infinite process spawning.
    main()