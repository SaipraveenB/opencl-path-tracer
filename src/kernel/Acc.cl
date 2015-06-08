//Implementing k-D Trees

typdef struct _box box;
struct box
{
    float2 xBounds;
    float2 yBounds;
    float2 zBounds;

    float4 centre;
};

typdef struct _node node;
struct node
{
    float splitVal;                 //Value used to split the volume
    int axis;                       //Axis used to split the volume
    node* left;                     //left subtree
    node* right;                    //right subtree
    int index;                      //index of the object (-1 if not leaf)
};

node newNode(float splitVal, int axis)
{
    node Node;
    Node.splitVal = splitVal;
    Node.axis = axis;
    return Node;    
}

node newLeaf(int index)
{
    node new;
    new.left = NULL;
    new.right = NULL;
    new.index = index;
    return new;
}

typedef struct _kdtree tree;
struct tree
{
    node* root;
    int count;
};

//TODO:: Figure out logic for intersection.
bool intersect();

//TODO:: Figure out logic for building the tree.
void buildTree();

//TODO:: Figure out logic to sort th objects.
