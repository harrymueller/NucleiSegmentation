using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class Generation : MonoBehaviour
{
    // local variables
    private Vector3 minPos, maxPos, sizes; 
    private float[] scale;
    private int numSpheres;
    private Material mat;
    private GameObject parent;
    private DensityDistribution dist;

    //private Color[] all_colours;
    private uint colour_i = 0;

    /*
    *   Setups positions, creates parent game object, and saves the required parameters to local variables
    */
    public void create(Vector3 sizes, int borderSize, float[] scale, int numSpheres, Material mat, DensityDistribution dist) {
        Vector3 sphereSizes = new Vector3(scale[1], scale[1], scale[1]);
        Vector3 borderSizes = new Vector3(borderSize, 0, borderSize);

        minPos = (sizes - sphereSizes) / -2.0f + borderSizes;
        maxPos = (sizes - sphereSizes) /  2.0f - borderSizes;

        this.sizes = sizes;
        this.scale = scale;
        this.numSpheres = numSpheres;
        this.mat = mat;

         // find a, b, c coefficients
        this.dist = dist;
        this.dist.calculate(minPos.x, maxPos.x);

        parent = new GameObject("Nuclei");
    }

    /*
    * Randomly creates nuclei
    */
    public void populate()
    {
        // delete nuclei if any exist already
        if (parent.transform.childCount > 0) empty();
        int i = 0;
        int failed = 0; // safety net

        while (i < numSpheres && failed < numSpheres) {
            if (createObj()) i++;
            else failed++;
        }
    }

    /*
    *   Deletes nuclei in the scene
    */
    public void empty() {
        GameObject.Destroy(parent);
        parent = new GameObject("Nuclei");
        colour_i = 0;
    }

    /*
    *   Checks whether a new object will have overlap with other objects
    */
    private bool CheckOverlap(GameObject obj, PrimitiveType type)
    {
        obj.SetActive(false);
        if (type == PrimitiveType.Sphere) { // sphere
            obj.GetComponent<SphereCollider>().isTrigger = true;
            return Physics.CheckSphere(obj.transform.position, obj.GetComponent<SphereCollider>().radius); 
        } else { // capsule
            // get start and end coords
            obj.GetComponent<CapsuleCollider>().isTrigger = false;
            Vector3 start; Vector3 end; float rad; 
            // https://github.com/justonia/UnityExtensions/blob/master/PhysicsExtensions.cs
            PhysicsExtensions.ToWorldSpaceCapsule(obj.GetComponent<CapsuleCollider>(), out start, out end, out rad);
            return Physics.CheckCapsule(start, end, rad);
        }
    }

    /*
    *   Randomly returns either a sphere or capsule primitive type
    */
    private PrimitiveType getRandomObj() {
        switch (Random.Range(0, 2)) {
            case 0:
                return PrimitiveType.Sphere;
            case 1:
                return PrimitiveType.Capsule;
            default:
                return PrimitiveType.Sphere;
        }
    }

    /*
    *   Generates a random Vector3 based on passed values
    */
    private Vector3 randomGenVector(float minVal, float maxVal) {
        Vector3 rand = new Vector3(0,0,0);
        for (int i = 0; i < 3; i++) rand[i] = Random.Range(minVal, maxVal);
        return rand;
    }

    /*
    *   Creates a random object
    */
    private bool createObj() {
        PrimitiveType objType = getRandomObj();
        GameObject obj = GameObject.CreatePrimitive(objType);

        // parent and material
        obj.transform.parent = parent.transform;
        obj.GetComponent<MeshRenderer>().material = mat;

        // add position, scale, rotation
        obj.transform.position = new Vector3(dist.convert(Random.Range(minPos.x, maxPos.x)), Random.Range(minPos.y, maxPos.y), Random.Range(minPos.z, maxPos.z));
        obj.transform.localScale = randomGenVector(scale[0], scale[1]);
        obj.transform.rotation = Quaternion.Euler(randomGenVector(0f, 360f));

        // shading
        float colour = (sizes.y + obj.transform.position.y) / sizes.y;
        obj.GetComponent<MeshRenderer>().material.SetColor("_Color", new Color(colour, colour, colour, 1f));

        // if obj not a plane and has no overlap, return true
        if (objType != PrimitiveType.Plane && CheckOverlap(obj, objType)) {
            GameObject.Destroy(obj);
            return false;
        } else {
            obj.SetActive(true);
            return true;
        }
    }

    /*
    *   For each nuclei, change the colour to a random one
    */
    public void randomiseColours() {
        Color c;
        foreach (Transform child in parent.transform) {
            // based on the color index, uses sets of 7bits to create a colour from (0.5, 1.0]
            float r = ((float) (colour_i>>14)      + 128) / 255; 
            float g = ((float) (colour_i>>7) % 128 + 128) / 255; 
            float b = ((float) (colour_i)    % 128 + 128) / 255; 

            c = new Color(r,g,b);
            child.gameObject.GetComponent<MeshRenderer>().material.SetColor("_Color", c);
            colour_i += 2;
        }
    }
}
