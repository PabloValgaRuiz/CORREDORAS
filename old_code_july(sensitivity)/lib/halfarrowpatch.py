import numpy as np
from matplotlib.patches import FancyArrowPatch

# Import directly from iplotx
from iplotx.edge.arrow import make_arrow_patch


class HalfArrowPatch(FancyArrowPatch):
    """
    A patch that draws a half arrow using iplotx's arrow implementation.
    
    This properly uses iplotx as a dependency rather than copying their code.
    """
    
    def __init__(self, posA, posB, marker="|/", width=8, height=None, curvature=None, **kwargs):
        """
        Create a half arrow from posA to posB using iplotx's arrow styles.
        
        Parameters
        ----------
        posA : tuple
            Starting position (x, y)
        posB : tuple
            Ending position (x, y) - arrow tip will be here
        marker : str
            Arrow marker style from iplotx. Options include:
            "|>", "|/", "|\\", ">", "<", ">>", ")>", ")", "(", "]", "[", "|",
            "x", "s", "d", "p", "q"
        width : float
            Width of the arrow in points
        height : float or None
            Length of the arrowhead in points. If None, uses width * 1.3
        **kwargs
            Additional arguments for line styling (color, linewidth, linestyle, alpha, etc.)
        """
        self.arrow_marker = marker
        self.arrow_width = width
        self.arrow_height = height
        
        # Extract arrow-specific kwargs
        self.arrow_color = kwargs.pop('arrow_color', kwargs.get('color', None))
        self.arrow_linewidth = kwargs.pop('arrow_linewidth', kwargs.get('linewidth', 1))
        self.arrow_alpha = kwargs.pop('arrow_alpha', kwargs.get('alpha', 1.0))
        
        # Don't use FancyArrowPatch's arrowstyle
        kwargs.pop('arrowstyle', None)
        
        # Handle curvature (similar to arc3)
        if curvature is not None:
            kwargs['connectionstyle'] = f'arc3,rad={curvature}'

        super().__init__(posA, posB, arrowstyle='-', **kwargs)

        # super().__init__(posA, posB, arrowstyle='-', **kwargs)
        
    def _get_iplotx_arrow_patch(self):
        """Get the arrow patch from iplotx."""
        kwargs = {}
        
        if self.arrow_color is not None:
            kwargs['color'] = self.arrow_color
        if self.arrow_linewidth is not None:
            kwargs['linewidth'] = self.arrow_linewidth
        if self.arrow_alpha is not None:
            kwargs['alpha'] = self.arrow_alpha
        
        # Use iplotx's make_arrow_patch function
        patch, size_max = make_arrow_patch(
            marker=self.arrow_marker,
            width=self.arrow_width,
            height=self.arrow_height,
            **kwargs
        )
        
        return patch
    
    def draw(self, renderer):
        """Override draw to add iplotx arrow at the end."""
        if not self.get_visible():
            return

        # Draw the main connecting line (with possible curvature)
        super().draw(renderer)

        if self._posA_posB[0] is None:
            return

        path = self.get_path()
        if path is None or len(path.vertices) < 2:
            return

        vertices = path.vertices
        tip = vertices[-1]
        prev = vertices[-2]

        dx = tip[0] - prev[0]
        dy = tip[1] - prev[1]
        angle = np.degrees(np.arctan2(dy, dx))

        # Get iplotx arrow patch
        arrow_patch = self._get_iplotx_arrow_patch()

        import matplotlib.transforms as mtransforms

        # Convert arrow size from points to pixels
        dpi = renderer.points_to_pixels(1.0)
        scale = dpi / 72.0

        # Compute offset so the arrow tip lands at posB ---
        path_bbox = arrow_patch.get_path().get_extents()
        tip_offset = path_bbox.xmax  # arrow tip in marker coordinates
        # Move origin so the tip sits at (0,0)
        trans = (mtransforms.Affine2D()
                .translate(-tip_offset, 0)  # shift left so tip is origin
                .scale(scale)
                .rotate_deg(angle)
                .translate(*tip))
        trans += self.get_transform()

        arrow_patch.set_transform(trans)
        arrow_patch.set_zorder(self.get_zorder())
        arrow_patch.draw(renderer)