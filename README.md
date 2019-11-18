## File Hierarchy
### Cluster specific configurations 
Cluster specific Makefile and slurm scripts are store in `./cluster_spec` directory.
#### For adroit
OSMesa library is installed and does not need to specify link.
#### For tiger
~~OSMesa library is located at ~/visit2_13_3.linux-x86_64~~
OSMesa library is installed and does not need to specify link.
#### For Della
OSMesa library is installed and does not need to specify link. In fact they are installed under the same directory `/lib86`
#### For stampede2
OSMesa library is installed under `$WORK/software` and needs to be linked.

### Project specific source code
project specific source code should be synced with desktop regularly.