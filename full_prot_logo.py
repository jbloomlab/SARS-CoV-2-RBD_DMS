import collections
import itertools
import math

import logomaker

import matplotlib.pyplot as plt

import pandas as pd


def full_prot_logo(
        tidy_df,
        *,
        site_col='site',
        letter_col='letter',
        height_col='height',
        color_col='color',
        highlight_color_col=None,
        highlight_alpha_col=None,
        sites_per_line=100,
        scalewidth=1,
        scaleheight=1,
        fade_letters_by_height=None,
        logo_kwargs=None,
        ylims=None,
        all_letters=('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'),
        missing_letter='error',
        letters_to_drop=('*',),
        style_xticks_kwargs=None,
        xlabel=None,
        ylabel=None,
        label_fontsize=16,
        baseline_on_top=True,
        ):
    """Draw logo wrapping several lines with custom colors and overlays.
    
    Parameters
    -----------
    tidy_df : pandas.DataFrame
        Holds data in tidy format, one line per letter.
    site_col : str
        Column in `tidy_df` with site number.
    letter_col : str
        Column in `tidy_df` with letter identity (e.g., amino acid).
    height_col : str
        Column in `tidy_df` with letter height.
    color_col : str
        Column in `tidy_df` with letter color.
    highlight_color_col : None or str
        Column in `tidy_df` with background highlight color, or `None`
        or `NA` if site not highlighted. Only one color can be assigned
        per site.
    highlight_alpha_col : None or float
        Column in `tidy_df` with background color alpha (transparency) for
        site highlighting. If not present, defaults to 0.25.
    sites_per_line : int
        Number of sites per line.
    scalewidth : float
        Scale overall figure height by this much.
    scaleheight : float
        Scale overall figure width by this much.
    fade_letters_by_height : None or 2-tuple
        If not `None`, set alpha transparency of letters proportional to
        their height, going from `(min_alpha, max_alpha)`.
    logo_kwargs : None or dict
        Keyword arguments to ``logomaker.Logo``. Key ones include 'width',
        'vpad', and 'font_name'.
    ylims : 2-tuple or None
        Y-axis limits, or `None` to auto-determine.
    all_letters : tuple or list
        All letters for which we plot heights.
    missing_letter : {'zero_height', 'error'}
        If letter is missing at a site, assign zero height or raise error?
    letters_to_drop : tuple or list
        Do not plot these letters.
    style_xticks_kwargs : None or dict
        Keyword arguments to pass to ``logomaker.Logo.style_xticks``. For
        instance, to change spacing between tick labels to every 10th site, use
        ``style_xticks_args={'spacing': 10}``.
    xlabel : str or None
        Label for x-axis (shared over entire plot).
    ylabel : str or None
        Label for y-axis (shared over entire plot).
    label_fontsize : int
        Size of labels drawn for `xlabel`, `ylabel`.
    baseline_on_top : bool
        Draw baseline (horizontal line at 0 height) on top of letters.

    """
    expect_cols = [site_col, letter_col, height_col, color_col]
    for col in [highlight_color_col, highlight_alpha_col]:
        if col is not None:
            expect_cols.append(col)
    for col in expect_cols:
        if col not in tidy_df.columns:
            raise ValueError(f"`tidy_df` lacks column {col}")
            
    if set(letters_to_drop).intersection(set(all_letters)):
        raise ValueError('overlap between `letters_to_drop` and `all_letters`')

    # make wide data frame for logomaker
    if tidy_df[site_col].dtype != int:
        raise ValueError('`site_col` currently must be integers')
    tidy_df = (tidy_df
               .query(f"{letter_col} not in {letters_to_drop}")
               .sort_values(site_col)
               )
    extra_letters = set(tidy_df[letter_col].unique()) - set(all_letters)
    if extra_letters:
        raise ValueError(f"`tidy_df` has extra letters: {extra_letters}\n"
                         'add them to `all_letters` or `letters_to_drop`')
    dup_entries = (
        tidy_df
        .groupby([site_col, letter_col])
        .aggregate(n_entries=pd.NamedAgg(height_col, 'count'))
        .query('n_entries > 1')
        )
    if len(dup_entries):
        raise ValueError(f"duplicates for sites / letters:\n{dup_entries}")
    wide_df = tidy_df.pivot_table(
                        index=site_col,
                        columns=letter_col,
                        values=height_col,
                        )
    if wide_df.isnull().any(None):
        if missing_letter == 'zero_height':
            wide_df = wide_df.fillna(0)
        elif missing_letter == 'error':
            raise ValueError('some letters missing, set `missing_letter`')
        else:
            raise ValueError(f"invalid `missing_letter` of {missing_letter}")

    # dict matching (site, letter) to color
    colors = tidy_df.set_index([site_col, letter_col])[color_col].to_dict()

    # dicts matching site to highlight color and alpha
    highlight_colors = {}
    highlight_alphas = collections.defaultdict(lambda: 0.25)
    for d, col in [(highlight_colors, highlight_color_col),
                   (highlight_alphas, highlight_alpha_col),
                   ]:
        if col is not None:
            site_vals = (tidy_df
                         [tidy_df[col].notnull()]
                         [[site_col, col]]
                         .drop_duplicates()
                         )
            dup_site_vals = (site_vals
                             .groupby(site_col)
                             .aggregate(n=pd.NamedAgg(col, 'count'))
                             .query('n > 1')
                             )
            if len(dup_site_vals):
                raise ValueError(f"multiple {col} for sites:\n{dup_site_vals}")
            for k, v in site_vals.set_index(site_col)[col].to_dict().items():
                d[k] = v

    # set up figure
    nsites = len(wide_df)
    nlines = math.ceil(nsites / sites_per_line)
    sites_per_line = min(sites_per_line, nsites)  # reduce if needed
    fig = plt.figure(figsize=(scalewidth * sites_per_line * 0.3,
                              scaleheight * nlines * 1.75),
                     )

    # map letters to fading
    if fade_letters_by_height:
        letter_fading = {}
        min_alpha, max_alpha = fade_letters_by_height
        if not 0 <= min_alpha < max_alpha <= 1:
            raise ValueError('fade_letters_by_height must span non-zero'
                             'range betweeen 0 and 1')
        min_height = wide_df.abs().min().min()
        max_height = wide_df.abs().max().max()
        for site, letter in itertools.product(wide_df.index, wide_df.columns):
            abs_height = abs(wide_df.at[site, letter])
            norm_fade = (abs_height - min_height) / (max_height - min_height)
            assert 0 <= norm_fade <= 1, norm_fade
            fade = norm_fade * (max_alpha - min_alpha) + min_alpha
            assert min_alpha <= fade <= max_alpha
            letter_fading[(site, letter)] = fade

    # auto-determine y-axis limits
    ypad = 1.02
    if ylims is None:
        if all(tidy_df[height_col] >= 0):
            ymin = 0
            ymax = ypad * tidy_df[height_col].max()
        elif all(tidy_df[height_col] <= 0):
            ymax = 0
            ymin = ypad * tidy_df[height_col].min()
        else:
            ymax = tidy_df[height_col].max()
            ymin = tidy_df[height_col].min()
            ymax += ypad * (ymax - ymin)
            ymin -= ypad * (ymax - ymin)

    # arguments for xtick styling
    xticks_kwargs = {'spacing': 5,  # number every five sites
                     'rotation': 90,  # rotated tick marks
                     'fontdict': {'verticalalignment': 'top',
                                  'horizontalalignment': 'center',
                                  'fontsize': 10},
                     }
    if style_xticks_kwargs is not None:
        for key, val in style_xticks_kwargs.items():
            xticks_kwargs[key] = val

    # draw logos for each line of figure
    for iline in range(nlines):  # loop over lines
        df = wide_df.iloc[iline * sites_per_line: (iline + 1) * sites_per_line]
        isites = df.index.tolist()  # sites being plotted on this axis
        ax = plt.subplot2grid(shape=(nlines, sites_per_line),
                              loc=(iline, 0),
                              colspan=len(df),  # number of sites for this line
                              fig=fig,
                              )
        logo = logomaker.Logo(
                    df=df,
                    ax=ax,
                    **logo_kwargs,
                    )

        # color letters
        for site, letter in itertools.product(isites, all_letters):
            style_kwargs = {}
            if (site, letter) in colors:
                style_kwargs['color'] = colors[(site, letter)]
            if fade_letters_by_height:
                style_kwargs['alpha'] = letter_fading[(site, letter)]
            if style_kwargs:
                logo.style_single_glyph(p=site, c=letter, **style_kwargs)

        # highlight sites
        for site, highlight_color in highlight_colors.items():
            if site in isites:
                logo.highlight_position(p=site, 
                                        color=highlight_color,
                                        alpha=highlight_alphas[site],
                                        )

        # format axes and ticks
        logo.style_spines(visible=False)
        logo.style_xticks(**xticks_kwargs)
        ax.tick_params(axis='x',
                       length=0,  # no xtick lines
                       pad=0,  # no padding between xtick labels and axis
                       )
        ax.set_ylim(ylims)
        ax.set_yticks([])

        # draw baseline on top of letters?
        if baseline_on_top:
            logo.draw_baseline(zorder=1)

    # set figure-wide axis labels: https://stackoverflow.com/a/53172335
    if xlabel or ylabel:
        ax_fig = fig.add_subplot(111, facecolor='none', frameon=False)
        ax_fig.tick_params(labelcolor='none', top=False, bottom=False,
                           left=False, right=False)
        if xlabel:
            ax_fig.set_xlabel(xlabel, fontsize=label_fontsize, labelpad=10)
        if ylabel:
            ax_fig.set_ylabel(ylabel, fontsize=label_fontsize)

    fig.tight_layout(h_pad=1.5)

    return fig


if __name__ == '__main__':
    import doctest
    doctest.testmod()
