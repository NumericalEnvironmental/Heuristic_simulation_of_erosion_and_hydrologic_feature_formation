# Heuristic_simulation_of_erosion_and_hydrologic_feature_formation

![Preview](https://numericalenvironmental.files.wordpress.com/2019/01/hydro_system.jpg?w=1632)

This julia language script reads an initial N x N grid of elevations across an x-y plane and simulates the formation of erosional and hydrologic features in response to user-supplied input parameters. This is accomplished by delineating total drainage areas associated with every cell (based on matrix algebra), filling in grid depressions, and specifying local erosion rates that are functions of drainage area as well as land surface gradient. Additional information, and results for an example problem, are provided in my blog (https://numericalenvironmental.wordpress.com/2019/01/02/simulating-erosion-and-watershed-formation-on-a-synthetically-generated-surface/).

A second example problem, which includes a python script for additional post-processing of the erosion script output to model the formation of playa lakes, is presented in another blog post: https://numericalenvironmental.wordpress.com/2019/01/29/a-python-model-for-hydrologic-and-topographic-constraints-on-playa-lake-formation/.

The script requires julia 1.0 or above, including the DelimitedFiles and SparseArrays packages. Two text input files are also required: (1) a comma delimited initial x-y-z surface file, labeled “surface_input.csv”, and (2) a parameters file (“params.text”) to adjust performance characteristics - see internal documentation for specific details.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

