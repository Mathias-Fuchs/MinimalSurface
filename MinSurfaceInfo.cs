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
                // return Properties.Resources.illu4;
                return new System.Drawing.Bitmap("Resources/illu4.bmp");
                }
        }
        public override string Description
        {
            get
            {
                 return "This library computes minimal surfaces bounded by one or two closed input curves," +
                    "without resorting to mesh relaxation.";
            }
        }
        public override Guid Id
        {
            get
            {
                return new Guid("638FBE32-8B19-41CD-B786-F8CB7A20D55E");
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
                return "mathias@mathiasfuchs.de";
            }
        }
    }
}
