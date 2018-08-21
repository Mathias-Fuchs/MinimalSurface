using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace MinSurface
{
    public class MinSurfaceInfo : GH_AssemblyInfo
    {
        public override string Name
        {
            get
            {
                return "Mini";
            }
        }
        public override Bitmap Icon
        {
            get
            {
                return new System.Drawing.Bitmap("Resources/illu4.bmp");
                }
        }
        public override string Description
        {
            get
            {
                 return "This library computes minimal surfaces bounded by a closed input curve, without resorting to mesh relaxation.";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("832052de-9224-4ef9-aba0-dad6c996f7eb");
            }
        }
        public override string AuthorName
        {
            get
            {                
                return "Mathias Fuchs";
            }
        }
        public override string AuthorContact
        {
            get
            {
                return "mathias2975@gmail.com";
            }
        }
    }
}
