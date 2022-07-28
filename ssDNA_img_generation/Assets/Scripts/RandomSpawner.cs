using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using System;

/*
    - collisions <- fix
    - save true segments
    - comments for classes
    - github
*/

public class RandomSpawner : MonoBehaviour
{
    public Camera cam;

    public Material black;
    public Material white;

    public Vector3 sizes = new Vector3(15, 3, 10);
    public int borderSize = 5;
    public int numSpheres = 10;
    public float[] scale = {0.7f, 1.0f};

    public bool useStepWiseEquations = true;
    public float[] proportions = {0.25f, 0.5f}; // 25% of nuclei in 50% of image -> 75% in 50%
    public int numPhotos = 10;
    private int numTaken; // counter for number of images generated

    public string filepath;
    public string fileformat = "{0}/{1}_{2}.png"; // dir/nuclei_0

    private Generation generator;

    // Start is called before the first frame update
    void Start()
    {
        // setup distribution
        DensityDistribution dist = gameObject.AddComponent<DensityDistribution>();
        if (useStepWiseEquations) dist.set_stepwise(proportions);

        generator = gameObject.AddComponent<Generation>();
        generator.create(sizes, borderSize, scale, numSpheres, white, dist);
        
        // populate for the first time
        generator.populate();

        // create plane and fix screen aspect ratio
        Misc.createBackgroundPlane(sizes, black);
        Misc.fixScreenSize(sizes, cam);
        
        Camera.onPreRender += preRenderFunc;
        Camera.onPostRender += postRenderFunc;
    }

    /*
    *   Prerender callback func to empty then populate scene
    */
    void preRenderFunc(Camera c) {
        // checks if enough photos has been taken
        if (numTaken >= numPhotos*2) OnDestroy();
        else if (numTaken == 0) return;// start occurs after prerender
        
        // otherwise generator new positions
        else if (numTaken % 2 == 0) { 
            generator.empty();
            generator.populate();
        } else {
            generator.randomiseColours();
        }
    }

    /*
    * Take a screenshot at the appropiate dimensions after rendering
    */
    void postRenderFunc(Camera c) {
        int width = RenderTexture.active.width;
        int height = RenderTexture.active.height;

        var t = new Texture2D(width, height, TextureFormat.RGB24, true);
        t.ReadPixels(new Rect(0, 0, width, height), 0, 0);
        t.Apply();
        var bytes = t.EncodeToPNG();

        Destroy(t);
        string type = (numTaken % 2 == 0) ? "nuclei" : "truth";
        File.WriteAllBytes(string.Format(fileformat, filepath, type, (numTaken++)/2), bytes);
    }

    /*
    *   Removes camera callback functions
    */
    void OnDestroy()
    {
        Camera.onPreRender -= preRenderFunc;
        Camera.onPostRender -= postRenderFunc;
    }
}
