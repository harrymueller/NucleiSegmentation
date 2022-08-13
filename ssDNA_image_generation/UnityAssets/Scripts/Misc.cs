using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class Misc
{
    /**
    * creates a black plane at y=-10 with correct size for the background to the image
    */
    public static void createBackgroundPlane(Vector3 sizes, Material mat) {
        GameObject obj = GameObject.CreatePrimitive(PrimitiveType.Plane);
        obj.transform.position = new Vector3(0,-10,0);
        obj.transform.localScale = new Vector3(sizes.x/10f, 1, sizes.z/10f);
        obj.GetComponent<MeshRenderer>().material = mat;
    }

    /*
    *   Forces the screen to be of the correct dimension
    */
    public static void fixScreenSize(Vector3 sizes, Camera cam) {
        // https://gamedesigntheory.blogspot.com/2010/09/controlling-aspect-ratio-in-unity.html
        float targetaspect = sizes.x / sizes.z;

        // determine the game window's current aspect ratio
        float windowaspect = (float)Screen.width / (float)Screen.height;

        // current viewport height should be scaled by this amount
        float scaleheight = windowaspect / targetaspect;

        // if scaled height is less than current height, add letterbox
        if (scaleheight < 1.0f)
        {  
            Rect rect = cam.rect;

            rect.width = 1.0f;
            rect.height = scaleheight;
            rect.x = 0;
            rect.y = (1.0f - scaleheight) / 2.0f;
            
            cam.rect = rect;
        }
        else // add pillarbox
        {
            float scalewidth = 1.0f / scaleheight;

            Rect rect = cam.rect;

            rect.width = scalewidth;
            rect.height = 1.0f;
            rect.x = (1.0f - scalewidth) / 2.0f;
            rect.y = 0;

            cam.rect = rect;
        }
        
        cam.orthographicSize = sizes.z / 2f;
    }

    public static void setupOutputs(string dir, string[] foldernames) {
        checkEmpty(dir + "/" + foldernames[0]);
        checkEmpty(dir + "/" + foldernames[1]);
    }

    private static void checkEmpty(string dir) {
        // if dir exists, delete
        if(Directory.Exists(dir)) Directory.Delete(dir, true); 
        
        Directory.CreateDirectory(dir);
    }
}
