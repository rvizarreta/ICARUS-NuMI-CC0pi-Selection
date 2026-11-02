import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

def mark_pot(ax, exposure, horizontal=False, vadj=0) -> None:
    """
    Add the POT information to the plot. The POT information will
    either be added to the top right corner of the plot along the
    horizontal axis or the top left corner of the plot along the
    vertical axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to add the POT information to.
    exposure : float
        The POT value to use in the marking.
    horizontal : bool, optional
        A flag to indicate if the POT information should be added
        along the horizontal axis. The default is True.
    vadj : float, optional
        The vertical adjustment to use when adding the POT information.
        The default is 0.

    Returns
    -------
    None.
    """
    mag = int(np.floor(np.log10(exposure)))
    usepot = exposure/10**mag
    s = f'NuMI {usepot:.2f}'+f'$\\times 10^{{{mag}}}$ POT'
    xrange = ax.get_xlim()
    yrange = ax.get_ylim()
    if horizontal:
        usey = yrange[1] + 0.025*(yrange[1] - yrange[0]) + vadj*(yrange[1] - yrange[0])
        usex = xrange[1] - 0.02*(xrange[1] - xrange[0])
        ax.text(x=usex, y=usey, s=s, fontsize=10, color='black', horizontalalignment='right')
    else:
        usey = yrange[1] + 0.02*(yrange[1] - yrange[0])
        usex = xrange[0] - 0.12*(xrange[1] - xrange[0])
        ax.text(x=usex, y=usey, s=s, fontsize=10, color='black', verticalalignment='top', rotation=90)

def mark_containment(ax, text, vadj=0, inside_y=None, up_pts=0):
    """
    Add the muon-containment selection, in a box, at the top right of the
    plot, in the place where mark_pot puts the POT on the 1D histograms
    (same right edge, same baseline). Expects the bold mathtext label from
    Variable.containment_label.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to add the text to.
    text : str
        The text to draw.
    vadj : float, optional
        The vertical adjustment, as a fraction of the y range. The default
        is 0.
    inside_y : float, optional
        If given, the box is drawn inside the axes instead, with its top edge
        at this height (axes fraction) and its right edge on the right edge
        of the axes. Used where the top-right corner above the axes is already
        taken by the POT label. The default is None.
    up_pts : float, optional
        Extra upward offset in points (outside-the-axes placement only). Used
        to put the box on the first title line, above the POT label. The
        default is 0.

    Returns
    -------
    matplotlib.text.Annotation
        The text artist (its bbox patch is the box).
    """
    yrange = ax.get_ylim()
    usey = yrange[1] + 0.025*(yrange[1] - yrange[0]) + vadj*(yrange[1] - yrange[0])
    fontsize, pad, lw = 8, 0.35, 0.8
    # The right edge of the box is put on the right edge of the axes: x is
    # 1.0 in axes coordinates, shifted left by the box padding
    # (pad * fontsize points) and half the edge line width.
    shift = pad * fontsize + lw / 2
    if inside_y is not None:
        shift += 6   # box inside the axes: leave 6 pt between its edge and the axes frame
    # Same format as the containment box on the 1D histograms: 8 pt bold
    # text in a thin black box with a white fill.
    if inside_y is not None:
        xy, xycoords, va = (1.0, inside_y), 'axes fraction', 'top'
    else:
        xy, xycoords, va = (1.0, usey), ('axes fraction', 'data'), 'baseline'
    return ax.annotate(text, xy=xy, xycoords=xycoords,
                xytext=(-shift, 0 if inside_y is not None else up_pts), textcoords='offset points',
                fontsize=fontsize, color='black',
                horizontalalignment='right', verticalalignment=va,
                annotation_clip=False,
                bbox=dict(boxstyle=f'square,pad={pad}', fc='white', ec='black', lw=lw))

def mark_preliminary(ax, label, vadj=0, hadj=0, fontsize=10, va='bottom',
                     black_prefix=False, subtitle_fontsize=10) -> None:
    """
    Add a preliminary label to the plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to add the preliminary label to.
    label : str
        The label to add to the plot to indicate that the plot is
        preliminary.
    vadj : float, optional
        The vertical adjustment to use when adding the label. The
        default is 0.
    hadj : float, optional
        The horizontal adjustment to use when adding the label. The
        default is 0.
    fontsize : float, optional
        Font size of the label. The default is 10.
    va : str, optional
        Vertical alignment of the label. The default is 'bottom'.
    black_prefix : bool, optional
        If True and the label starts with a bold mathtext word (e.g.
        '$\\bf{ICARUS}$ Data'), that word is drawn in black and the rest of
        the label in the Data/MC color (both in the style sheet font). Both pieces share the
        same baseline. The default is False.
    subtitle_fontsize : float, optional
        Font size of the second line of the label when black_prefix is True.
        The default is 10 (the size of the POT label).

    Returns
    -------
    None.
    """
    yrange = ax.get_ylim()
    usey = yrange[1] + 0.025*(yrange[1] - yrange[0]) + vadj*(yrange[1] - yrange[0])
    xrange = ax.get_xlim()
    usex = xrange[0] + 0.025*(xrange[1] - xrange[0]) + hadj*(xrange[1] - xrange[0])
    color = 'forestgreen' if 'data' in label.lower() else 'blue'
    m = None
    if black_prefix:
        first, _, subtitle = label.partition('\n')
        m = re.match(r'^(\$\\bf\{([^}]*)\}\$)\s+(\S+)(?:\s+(.*))?$', first)
    if m is not None:
        try:
            fig = ax.figure
            renderer = fig.canvas.get_renderer()
            gap = 0.3 * fontsize
            # Two-line title. The (optional) second line sits on the usual
            # baseline, next to the POT label, and the main line is raised
            # above it by line_pts points.
            sub_fs = subtitle_fontsize
            # Baseline-to-baseline distance: cap height of the second line
            # (~0.72 em) + descender of the main line (~0.22 em) + 3 pt gap.
            # The Data/MC box hangs about 0.2 em below the main baseline, so
            # 2 more points are left between the lines.
            line_pts = (0.72 * sub_fs + 0.22 * fontsize + 5) if subtitle else 0
            if subtitle:
                ax.annotate(subtitle, xy=(usex, usey), xycoords='data',
                            fontsize=sub_fs, color='black',
                            horizontalalignment='left', verticalalignment=va,
                            annotation_clip=False)
            # Main line: "ICARUS" (bold, black), "Data"/"MC" (bold, colored),
            # then the rest (regular, black). Plain text, not mathtext, so
            # that it uses the style sheet's sans-serif font like the axis
            # labels. Each piece is anchored to the bottom-right corner of the
            # previous one (resolved at draw time, so the gaps are exact);
            # that corner lies one descent below the previous baseline, so
            # the descent is added back to keep all pieces on one baseline.
            head = ax.annotate(m.group(2), xy=(usex, usey), xycoords='data',
                               xytext=(0, line_pts), textcoords='offset points',
                               fontsize=fontsize, fontweight='bold', color='black',
                               horizontalalignment='left', verticalalignment=va,
                               annotation_clip=False)
            descent_px = renderer.get_text_width_height_descent(
                'lp', head.get_fontproperties(), False)[2]
            descent_pts = descent_px * 72.0 / renderer.dpi

            def _next(prev, text, extra_x=0, descent=None, size=None, **kw):
                return ax.annotate(text, xy=(1, 0), xycoords=prev,
                                   xytext=(gap + extra_x, descent_pts if descent is None else descent),
                                   textcoords='offset points',
                                   fontsize=fontsize if size is None else size,
                                   horizontalalignment='left', verticalalignment='baseline',
                                   annotation_clip=False, **kw)

            # Data/MC: white bold text on a rounded box filled with the
            # Data/MC color. The box extends pad_pts beyond the text on each
            # side, so the neighbouring gaps are widened by the same amount.
            # The word is drawn smaller than the rest of the line so that its
            # box fits the line height. Pieces are chained by bbox corners,
            # which lie one descent below their own baseline, so each step
            # adds back the descent of the piece it is anchored to.
            word_size = 0.8 * fontsize
            pad = 0.2
            pad_pts = pad * word_size
            # Vertical centering: the box (text bbox + pad) is centered at
            # (h - 2d)/2 above the word's baseline, where h and d are the line
            # height and descent of the word's font. The cap height of DejaVu
            # Sans is ~0.73 em, so the centre of the "ICARUS" capitals sits
            # 0.365*fontsize above the shared baseline. The word is raised by
            # the difference so that box centre and capital centre coincide.
            wfp = head.get_fontproperties().copy()
            wfp.set_size(word_size)
            _, w_h_px, w_d_px = renderer.get_text_width_height_descent('lp', wfp, False)
            box_center_pts = (w_h_px - 2 * w_d_px) / 2.0 * 72.0 / renderer.dpi
            shift = 0.365 * fontsize - box_center_pts
            word = _next(head, m.group(3), extra_x=pad_pts, size=word_size,
                         descent=descent_pts + shift,
                         color='white', fontweight='bold',
                         bbox=dict(boxstyle=f'round,pad={pad},rounding_size=0.5',
                                   fc=color, ec=color))
            if m.group(4):
                word_descent_px = renderer.get_text_width_height_descent(
                    'lp', word.get_fontproperties(), False)[2]
                _next(word, m.group(4), extra_x=pad_pts,
                      descent=word_descent_px * 72.0 / renderer.dpi - shift, color='black')
            return
        except Exception:
            pass
    ax.text(x=usex, y=usey, s=label, fontsize=fontsize, color=color, verticalalignment=va)

def draw_error_boxes(ax, x, y, xerr, yerr, **kwargs):
    """
    Adds error boxes to the input axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to which the error boxes are to be added.
    x : numpy.array
        The x-coordinates of the error boxes.
    y : numpy.array
        The y-coordinates of the error boxes.
    xerr : numpy.array
        The x-error values of the error boxes.
    yerr : numpy.array
        The y-error values of the error boxes.
    kwargs: dict
        Keyword arguments to be passed to the errorbar function.

    Returns
    -------
    None.
    """
    boxes = [Rectangle((x[i] - xerr[i], y[i] - yerr[i]), 2 * np.abs(xerr[i]), 2 * yerr[i]) for i in range(len(x))]
    pc = PatchCollection(boxes, **kwargs)
    ax.add_collection(pc)
