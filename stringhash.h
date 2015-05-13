
/* This stringhash is a simple tree structure where each node corresponds to a letter of a key string and a value is a nonnegative integer (-1 is used for errors, so the stored value must be nonnegative). It does not support individual key deletion. Over-writing an existing key is supported. */


typedef struct hashtreenode
{
   char letter;
   int val;
   int size;
   struct hashtreenode** next;
} treenode;

treenode* insert_key(char*,int,treenode*);
int search_treehash(char*,treenode*);
treenode* new_treenode(char);
treenode* new_lasttreenode(char,int);
treenode* new_tree(void);
treenode* new_branch(treenode*);
void delete_treenodeptrarray(treenode**, int);
void delete_treenode(treenode*);

