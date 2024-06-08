#' @name DuplexTrack
#' @title  class for the visualization of RNA duplexes

#' @description
#' Inherits the \code{Gviz::AnnotationTrack}, plots interaction ranges as boxes.
#' Arguments from \code{Gviz::AnnotationTrack}, as `stacking` which set boxes
#' layout are accepted. Parent aesthetics for labels are overwritten with Display
#' parameters of this class.
#' Accepts `GInteractions` object to plot and `GRanges` to define plot region
#'
#' Duplexes which can be displayed on the plot range are connected with arcs.
#' Duplexes which are partially outside of the range are displayed without arcs.
#' Labeles and appearance can be controlled with display parameters
#'
#' @param gi An `GInteractions` object
#' @param gr_region `GRanges` region for plotting
#' @param from Integer start coordinate of subset region. Used if `gr_region` is not provided
#' @param to Integer end coordinate of subset region. Used if `gr_region` is not provided
#' @param chromosome Chromosome of subset region. Used if `gr_region` is not provided
#' @param strand Used if `gr_region` is not provided
#' @param fill.column used for fill. Default is "" (empty) and triggers IGV color pallete.
#' **Display parameters**
#' \describe{
#'   \item{arcs.color}{Character. Color of the arcs. Default is "black".}
#'   \item{arc.location}{Character in c('inner','outer','midpoint'). Location of the arcs in X axis relative to range.
#'   Default is "inner"}
#'   \item{labels.v.offset.base}{Numeric. Base vertical offset for the labels. Default is 0.2. Other offesets are added to it.
#'   }
#'   \item{labels.v.offset.trans}{Numeric. Vertical offset for trans labels. Applied when one part of the duplex is outside of the plot.
#'    Recommended ranges are in -0.5 to 0.5 Default is 0.0.   }
#'   \item{labels.h.offset.trans}{Numeric. Horizontal offset for trans labels.
#'    Applied when one part of the duplex is outside of the plot
#'   Value is in nucleotide units. Default is 0.0. }
#'   \item{labels.v.offset.cis}{Numeric. Vertical offset for cis labels.  Recommended ranges are in -0.5 to 0.5 Default is 0.0. Default is 0.0. }
#'   \item{labels.h.offset.cis}{Numeric. Horizontal offset for cis labels.  Value is in nucleotide units. Default is 0.0.}
#'   \item{labels.fontsize}{Numeric. Font size of the labels. Default is 18.}
#'   \item{label.cis.above}{Logical. Whether the cis labels should be above. When set to FALSE, labels are plot for each box separately.
#'   Default is TRUE}
#'   \item{annotation.column1}{Character. First annotation column to use for labels. Default is "group" and generated internally.}
#'   \item{annotation.column2}{Character. Second annotation column to use for labels. Default is "" (empty).}
#'   \item{fill.column}{Character. Column used for fill. Default is "" (empty) and triggers IGV color pallete.}
#'   \item{labels.color}{Character. Color of the labels. Default is 'black'.}
#'   \item{labels.align}{Character. Alignment of the labels. Default is 'center'.
#'   Possible values are in c('left','right','center)}
#'   \item{arcConstrain}{Numeric. Minimum gap distance between arms of the
#'    interaction to draw arcs}
#' }
#' @export DuplexTrack
#' @examples
#' library(InteractionSet)
#' library(Gviz)
#' # generate input
#' anchor1 <- GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges(
#'         start = c(100, 600, 1100, 1600, 2100, 150, 400),
#'         end = c(200, 700, 1200, 1700, 2200, 250, 500)
#'     ),
#'     strand = "+"
#' )
#' anchor2 <- GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges(
#'         start = c(300, 800, 1300, 1800, 2300, 1500, 1700),
#'         end = c(400, 900, 1400, 1900, 2400, 1600, 1800)
#'     ),
#'     strand = "+"
#' )
#'
#' interactions <- GInteractions(anchor1, anchor2, mode = "strict")
#' # define plotting range
#' gr_region <- range(anchor1, anchor2)
#' interactions$anno_A <- sample(LETTERS, length(interactions))
#' interactions$anno_B <- interactions$anno_A
#' a <- DuplexTrack(interactions, gr_region = gr_region, stacking = "dense")
#' plotTracks(a, stacking = "dense")
#' plotTracks(a, stacking = "squish", annotation.column1 = "anno_A")
#'
#' # add interactions which are not fully in plot range: outside the range or on different chromosome()
#'
#' # one left (A) interaction arm outside of the plot, other on different chromosome
#' new_anchor1 <- GRanges(
#'     seqnames = c("chr1", "chr2"),
#'     ranges = IRanges(
#'         start = c(10, 600),
#'         end = c(90, 700)
#'     ),
#'     strand = "+"
#' )
#' new_anchor2 <- GRanges(
#'     seqnames = c("chr1", "chr1"),
#'     ranges = IRanges(
#'         start = c(1500, 1000),
#'         end = c(1600, 1200)
#'     ),
#'     strand = "+"
#' )
#'
#' new_interactions <- GInteractions(new_anchor1, new_anchor2)
#' new_interactions$anno_A <- c("A.out", "A.out_chr")
#' new_interactions$anno_B <- c("B.in", "B.in")
#' all_interactions <- c(interactions, new_interactions)
#'
#' b <- DuplexDiscoverer::DuplexTrack(all_interactions,
#'     gr_region = gr_region,
#'     annotation.column1 = "anno_A",
#'     annotation.column2 = "anno_B"
#' )
#'
#' plotTracks(b)
#'
#' # to customize plot, one can call, to see options
#' DuplexDiscoverer::availableDisplayPars(b)
#'
setClass(
    "DuplexTrack",
    contains = c("AnnotationTrack"),
    representation = representation(
        giobject = "StrictGInteractions",
        name = "character",
        gr_region = "GenomicRanges"
        #      fill = "vector"
        #      variables = "list",
        #      chromosome = "character",
        #      start = "numeric",
        #      end = "numeric",
        #      stacking = "character",
        #      initval = "character",
        #      range = "GenomicRanges",
    ),
    prototype = prototype(
        name = "DuplexTrack",
        dp = DisplayPars(
            arcs.color = "black",
            arc.location = "inner",
            arcConstrain = 4,
            labels.v.offset.base = 0.2,
            labels.v.offset.trans = 0.0,
            labels.h.offset.trans = 0.0,
            labels.v.offset.cis = 0.0,
            labels.h.offset.cis = 0.0,
            labels.fontsize = 18,
            label.cis.above = TRUE,
            annotation.column1 = "group",
            annotation.column2 = "",
            fill.column = "",
            labels.color = "black",
            labels.align = "center"
        )
    )
)

setMethod("initialize", "DuplexTrack", function(.Object, ...) {
    if (is.null(list(...)$range) && is.null(list(...)$genome) && is.null(list(...)$chromosome)) {
        return(.Object)
    }
    ## the display parameter defaults
    # Gviz:::.makeParMapping()
    # .Object <- Gviz:::.updatePars(.Object, "DuplexTrack")
    # .gv_parmapping()
    .Object <- .gv_updatepars(.Object, "DuplexTrack")
    range <- list(...)$range
    giobject <- list(...)$giobject
    .Object@giobject <- giobject
    # .Object@fill <- list(...)$fill
    .Object@gr_region <- list(...)$gr_region

    .Object <- callNextMethod()
    message("init")

    return(.Object)
})

DuplexTrack <- function(
        gi, start = NULL, end = NULL, gr_region = NULL,
        group, id, strand, chromosome, fill = NULL, fill.column = "",
        genome, stacking = "squish", name = "DuplexTrack", selectFun, importFunction,
        stream = FALSE, ...) {
    message("constructor")
    ## Some defaults
    if (is.null(gr_region)) {
        if (!any(c(is.null(start), is.null(strand), is.null(strand)))) {
            gr_region <- GRanges(as.character(chromosome), IRanges(start, end), strand = strand)
        }
    } else {
        message("Using provided Granges for plot")
    }

    # first select gi where at least single region is on the range
    gi_sbs <- subsetByOverlaps(gi, gr_region)
    if (length(gi_sbs) == 0) {
        errorCondition(message("No overlap found for the object ot plot and defined plotting range"))
    }

    # load fill colors or use pallete
    if (!is.null(fill.column) & (fill.column %in% colnames(mcols(gi_sbs)))) {
        # use colum with fill
        fill <- mcols(gi)[, fill.column]
    } else {
        # Use igv pallete
        mypal <- pal_igv("default", alpha = 0.9)(length(gi_sbs))
        message("Using IGV pallete for fill")
        mypal <- ggsci::pal_igv("default", alpha = 0.9)(length(gi_sbs))
        fill <- mypal
    }
    range_meta <- tibble(gap = pairdist(gi_sbs, type = "gap"), ) %>%
        mutate(
            in_range = ifelse(is.na(gap), 0, 1),
            arms_overlap = ifelse(in_range == 1 & gap <= 0, 1, 0)
        ) %>%
        data.frame()
    gi_sbs$gap <- NULL
    gi_sbs$arms_overlap <- NULL
    gi_sbs$in_range <- NULL
    mcols(gi_sbs) <- cbind(data.frame(mcols(gi_sbs)), range_meta)
    # select those boxes which are on the plotting range
    gi_boxes <- convert_gi_to_ranges(gi_sbs)
    range <- subsetByOverlaps(gi_boxes, gr_region)
    if (length(range) == 0) {
        stop("No overlaps found")
    }
    add_meta <- as_tibble(mcols(range)[c("id", "group")]) %>%
        mutate(arm = str_split_i(id, "\\.", 2)) %>%
        group_by(group) %>%
        mutate(n = n()) %>%
        ungroup() %>%
        mutate(
            full_in_plot = if_else(n == 2, 1, 0),
            pair_arm = if_else(arm == "B", "A", "B")
        ) %>%
        dplyr::select(arm, full_in_plot, pair_arm) %>%
        data.frame()

    mcols(range) <- cbind(mcols(range), add_meta)
    range$feature <- "unknown"
    range$density <- 1

    if (missing(chromosome) || is.null(chromosome)) {
        chromosome <- if (length(range) > 0) .chrName(as.character(seqnames(range)[1])) else "chrNA"
    }
    ## And finally the object instantiation, we have to distinguish between DetailsAnnotationTracks and normal ones
    # genome <- Gviz:::.getGenomeFromGRange(range, ifelse(is.null(genome), character(), genome[1]))
    genome <- rtracklayer::genome(range)
    return(new("DuplexTrack",
        giobject = (gi), start, end,
        chromosome = as.character(chromosome[1]), range = range,
        gr_region = gr_region,
        fill = fill, groupAnnotation = "group",
        name = name, genome = genome, stacking = stacking, ...
    ))
}


#' Show method for DuplexTrack
#'
#' @param object DuplexTrack.
#' @importFrom methods slotNames show
#' @return class representation
#' @export
#'
#' @examples
#' library(InteractionSet)
#' anchor1 <- GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges(
#'         start = c(100, 600, 1100, 1600, 2100),
#'         end = c(200, 700, 1200, 1700, 2200)
#'     ),
#'     strand = "+"
#' )
#' anchor2 <- GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges(
#'         start = c(300, 800, 1300, 1800, 2300),
#'         end = c(400, 900, 1400, 1900, 2400)
#'     ),
#'     strand = "+"
#' )
#'
#' interactions <- GInteractions(anchor1, anchor2, mode = "strict")
#' gr_region <- range(anchor1, anchor2)
#' a <- DuplexTrack(interactions, gr_region = gr_region, stacking = "dense")
#' show(a)
setMethod("show", signature(object = "DuplexTrack"), function(object) {
    slots <- methods::slotNames(object)
    cat(names(object), " class")
    selected_slots <- slots[slots != "dp"]
    for (slot_name in selected_slots) {
        cat((slot_name), ": \n")
        show(slot(object, slot_name))
    }
    pars_to_display <- names(availableDisplayPars(class(object)))
    plist <- object@dp@pars[pars_to_display]
    cat("Display parameters: \n")
    cat(paste(names(plist), plist, sep = ": ", collapse = "\n"))
})


#' Draw methods for DuplexTrack
#'
#' `Gviz::AnnotationTrack` stacking algorithm is used to calculate vertical
#' distribution of boxes for the interactions.
#' Boxes coordinates are later imported for placing labels and arcs
#' @keywords internal
#' @returns pushes boxes, arcs and labels to viewport
setMethod("drawGD", signature("DuplexTrack"), function(GdObject, minBase, maxBase, prepare = FALSE, subset = TRUE, ...) {
    imageMap(GdObject) <- NULL
    if (!length(GdObject)) {
        return(invisible(GdObject))
    }
    # Brought here from AnnotationTrack
    ## In prepare mode we need to make sure that the stacking information is updated from the optional display parameter (by calling
    ## the StackedTrack drawGD method) and also perform the collapsing of track items which could potentially lead to re-stacking.
    if (prepare) {
        message("prepare")
        pushViewport(viewport(xscale = c(minBase, maxBase), yscale = c(0, 1)))
        popViewport(1)
        return(invisible(GdObject))

        message("P. start: ", minBase)
        message("P. end: ", maxBase)
    }
    ## If there are too many stacks for the available device resolution we cast an error, otherwise we set the viewport
    bins <- stacks(GdObject)
    stacks <- max(bins)
    rev <- Gviz:::.dpOrDefault(GdObject, "reverseStrand", FALSE)

    temp_info <- GdObject@range

    # vector defining ids which are present in both reads
    vec_both_present <- as.integer(names(table(temp_info$group)[which(table(temp_info$group) == 2)]))
    # set the row of reads
    mcols(temp_info)$row_plot_id <- as.vector(bins)
    # provide info if gap is sufficient for arcs
    arcConstrain <- Gviz:::.dpOrDefault(GdObject, "arcConstrain")
    mcols(temp_info)$distance_gap <- ifelse(mcols(temp_info)$gap > arcConstrain, 1, 0)
    mcols(temp_info)$distance_gap[is.na(mcols(temp_info)$distance_gap)] <- 0
    # provide info about height required for blocks
    vec_double_height <- c()
    for (internal in seq_len(stacks)) {
        tmp <- temp_info[which(mcols(temp_info)$row_plot_id == internal)]
        unique_dist <- unique(mcols(tmp)$distance_gap)
        vec_double_height <- c(vec_double_height, if (length(unique_dist) == 1) unique_dist else 1)
    }
    vec_double_height_axis <- vec_double_height + 1
    vec_double_text <- 1 - vec_double_height

    # create boxes height
    # for usual box we will use 1 and for for boxes with requirements of arc 2
    # plus we should have text size
    annot_size <- Gviz:::.dpOrDefault(GdObject, "labels.v.offset.base")
    xscale <- if (!rev) c(minBase, maxBase) else c(maxBase, minBase)

    height_formula <- sum(vec_double_height_axis) + sum(vec_double_text * 0.05) + 1
    yscale <- if (!Gviz:::.dpOrDefault(GdObject, "reverseStacking", FALSE)) c(1, height_formula) else c(height_formula, 1)

    pushViewport(dataViewport(xscale = xscale, extension = 0, yscale = yscale, clip = TRUE))
    res <- Gviz:::.pxResolution(coord = "x")
    curVp <- Gviz:::vpLocation()

    # Q: I do not think we weed this warning
    if (curVp$size["height"] / stacks < Gviz:::.dpOrDefault(GdObject, "min.height", 3)) {
        stop("Too many stacks to draw. Either increase the device size or limit the drawing to a smaller region.")
    }

    ## We adjust the color saturation to indicate overplotting if necessary
    if (Gviz:::.dpOrDefault(GdObject, "showOverplotting", FALSE)) {
        dens <- as.numeric(values(GdObject)$density)
        if (length(unique(dens)) != 1) {
            minSat <- max(0.25, 1 / max(dens))
            minDens <- min(dens)
            saturation <- minSat + ((dens - minDens) / rDens / (1 / (1 - minSat)))
            bc <- unique(.getBiotypeColor(GdObject))
            baseCol <- rgb2hsv(col2rgb(bc))
            desatCols <- unlist(lapply(saturation, function(x) hsv(baseCol[1, ], x, baseCol[3, ])))
            names(desatCols) <- paste(unique(feature(GdObject)), rep(dens, each = length(bc)), sep = "_")
            feature(GdObject) <- paste(feature(GdObject), dens, sep = "_")
            desatCols <- desatCols[unique(names(desatCols))]
            displayPars(GdObject) <- as.list(desatCols)
        }
    }
    # Starting with boxes coordinates
    # bottom and top positions for every case
    internal <- 1
    vector_fin_bot <- c()
    vector_fin_top <- c()
    for (internal in seq(1, stacks)) {
        val_based_box_height <- seq(1, stacks)[internal]
        val_text_height <- sum((vec_double_text * annot_size)[seq_len(internal)])
        val_arc_height <- sum(vec_double_height[c(seq_len(internal))])
        val_fin_bot <- height_formula - val_based_box_height - val_text_height - val_arc_height
        val_fin_top <- val_fin_bot + 1
        vector_fin_bot <- c(vector_fin_bot, val_fin_bot)
        vector_fin_top <- c(vector_fin_top, val_fin_top)
        # print(paste(paste("Bottom coordinates row:",internal,"   "),height_formula,val_based_box_height,val_text_height,val_arc_height,val_fin_bot,sep="|"))
    }
    df_coord <- data.frame(row_plot_id = seq(1, stacks), top = vector_fin_top, bot = vector_fin_bot)

    ## Now we can pre-compute all the coordinates and settings for the elements to be drawn, ...
    box <- Gviz:::.boxes(GdObject, (stacks - bins) + 1)
    tmp_box <- merge(box, mcols(temp_info)[, c("row_plot_id", "id")], by = "id")
    tmp_box <- merge(tmp_box, df_coord, by = "row_plot_id")
    tmp_box$cy1 <- tmp_box$bot
    tmp_box$cy2 <- tmp_box$top
    tmp_box$textY <- tmp_box$top + annot_size
    row.names(tmp_box) <- tmp_box$id
    box <- tmp_box

    ## all the necessary display parameters from gviz
    barsAndLab <- Gviz:::.barsAndLabels(GdObject)
    bar <- barsAndLab$bars
    bartext <- barsAndLab$labels
    shape <- "box"

    # Take from AnnotationTrack
    col.line <- Gviz:::.dpOrDefault(GdObject, "col.line")[1]
    border <- Gviz:::.dpOrDefault(GdObject, "col")[1]
    lwd <- Gviz:::.dpOrDefault(GdObject, "lwd", 2)
    lty <- Gviz:::.dpOrDefault(GdObject, "lty", 1)
    alpha <- Gviz:::.dpOrDefault(GdObject, "alpha", 1)
    rotation <- Gviz:::.dpOrDefaultFont(GdObject, "rotation", "item", 0)
    rotation.group <- Gviz:::.dpOrDefaultFont(GdObject, "rotation", "group", 0)
    just <- Gviz:::.dpOrDefault(GdObject, "just.group", "above")

    # Our added dupslaypars
    color.arcs <- Gviz:::.dpOrDefault(GdObject, "arcs.color")
    font_gpar <- Gviz:::.fontGp(GdObject, subtype = "group")
    font_gpar$fontsize <- displayPars(GdObject, "labels.fontsize")
    font_gpar$col <- displayPars(GdObject, "labels.color")
    arc.location <- displayPars(GdObject, "arc.location")
    labels.v.offset.trans <- displayPars(GdObject, "labels.v.offset.trans")
    labels.h.offset.trans <- displayPars(GdObject, "labels.h.offset.trans")
    labels.v.offset.cis <- displayPars(GdObject, "labels.v.offset.cis")
    labels.h.offset.cis <- displayPars(GdObject, "labels.h.offset.cis")
    labels.align <- displayPars(GdObject, "labels.align")
    label.cis.above <- displayPars(GdObject, "label.cis.above")
    box.fill.colors <- displayPars(GdObject, "fill")
    # can be custom, but not tested for now
    height_factor_box <- 0.8
    height_factor_arc <- 0.8

    ## Plotting of the boxes
    message("Plotting Boxes")
    box$col <- if (is.na(border)) box$fill else border
    fill_vector <- box.fill.colors
    # assign correct fill colours
    if (!is.null(fill_vector)) {
        # Get unique indexes for each pair
        pair_indexes <- unique(sub("\\..*", "", rownames(box)))
        # Create a named vector with pairs and their corresponding colors
        pair_colors <- setNames(fill_vector, pair_indexes)
        # Map colors to each pair
        pair_colors_vector <- vapply(sub("\\..*", "", rownames(box)), function(index) pair_colors[[index]],"")
        box$fill <- as.vector(pair_colors_vector)
    }

    # adjust before plotting boxes
    box2 <- box
    hhe <- box2$cy2 - box2$cy1
    height_orig <- hhe[1]
    height_box <- hhe[1] * height_factor_box
    height_red <- height_orig - height_box
    box2$cy2 <- box2$cy1 + height_box
    box2$textY <- box2$textY - height_red

    # plot boxes
    .gv_plotboxes(box2, lwd = lwd, lty = lty, alpha = alpha)
    box <- box2

    message("Plotting Labels")
    # EGORS for compatibility tmp box is
    tmp_box <- box
    # create arcs
    message("Plotting arcs")
    # select unique ids with both arms available
    result <- as_tibble(tmp_box) %>%
        dplyr::filter(gap > arcConstrain) %>%
        group_by(group = sub("\\..*", "", id)) %>%
        summarize(has_A = any(grepl("A", id)), has_B = any(grepl("B", id))) %>%
        dplyr::filter(has_A & has_B)
    # if no ars are to plot - don't try to call graphics to avoid confusion
    if (nrow(result) == 0) {
        stop("Cannot plot provided region. Possible reason: interactions overlap.
           Suggestion: try to convert input to long GRanges and plot with
           standard Gviz Annotation track")
    }
    pass_arcs <- result$group

    # Egor S added df for max arc height
    arc_h_list <- list()
    df_arc_y <- data.frame()
    if (length(pass_arcs) != 0) {
        internal <- 1
        for (internal in seq_along(pass_arcs)) {
            grp_var <- paste(pass_arcs[internal], c("A", "B"), sep = ".")
            tmp <- tmp_box[grp_var, ]
            # options for creating arcs start end end coordinates
            # 1) midpoint
            # 2) inner
            # 3) outer
            if (arc.location == "midpoint") {
                points_x <- (abs(tmp$end) + abs(tmp$start)) / 2
                point_y <- max(tmp$cy2)
            } else if (arc.location == "inner") {
                points_x <- c(abs(tmp$end[1]), abs(tmp$start[2]))
                point_y <- max(tmp$cy2)
            } else if (arc.location == "outer") {
                points_x <- c(abs(tmp$start[1]), abs(tmp$end[2]))
                point_y <- max(tmp$cy2) + annot_size
            }
            # code for default
            y_base <- point_y #- height_red
            y_elev <- height_factor_arc * height_box

            r <- (points_x[2] - points_x[1]) / 2
            a <- y_elev / r / r
            x <- seq(-r, r, length.out = 20)
            y <- -a * x * x
            # shift
            x <- x + r + points_x[1]
            y <- y + y_elev + y_base
            arc_h_list[["grpvar"]] <- grp_var
            arc_h_list[["arc_max_y"]] <- max(y)
            arc_h_list[["arc_center_x"]] <- points_x[1] + (points_x[2] - points_x[1]) / 2
            df_arc_y <- rbind(df_arc_y, data.frame(arc_h_list))
            grid.xspline(x, y,
                shape = -1, open = TRUE,
                default.units = "native",
                gp = gpar(lwd = lwd, lty = lty, col = color.arcs)
            )
            ############################################
        }
        df_arc_y <- df_arc_y %>% mutate("grp" = str_split_i(grpvar, "\\.", 1))
    }
    # text exists and asked to plot
    # depending on the annotation annotation.column2, set labels for cis

    if (nrow(box) > 0) {
        tmp_box2 <- as_tibble(tmp_box) %>% mutate(coordinate = (start + end) / 2)
        if (Gviz:::.dpOrDefault(GdObject, "annotation.column2") == "") {
            annotation.column1 <- Gviz:::.dpOrDefault(GdObject, "annotation.column1")
            annotation.column2 <- Gviz:::.dpOrDefault(GdObject, "annotation.column1")
        } else {
            annotation.column1 <- Gviz:::.dpOrDefault(GdObject, "annotation.column1")
            annotation.column2 <- Gviz:::.dpOrDefault(GdObject, "annotation.column2")
        }
        tmp_box2 <- tmp_box2[c("id", "textY", "coordinate", annotation.column1, annotation.column2)]
        colnames(tmp_box2) <- c("id", "y", "x", "txt1", "txt2")
        # apply offest
        tmp_box2$yOffset.A <- tmp_box2$y
        tmp_box2$yOffset.B <- tmp_box2$y
        tmp_box2$arm <- str_split_i(tmp_box$id, "\\.", 2)
        tmp.A <- tmp_box2[str_detect(tmp_box2$id, ".A"), ]
        tmp.B <- tmp_box2[str_detect(tmp_box2$id, ".B"), ]
        # for per-box
        yOffset.A <- tmp.A$y
        yOffset.B <- tmp.B$y
        if (label.cis.above == TRUE) {
            # CIS LABELS
            if (nrow(df_arc_y) != 0) {
                cislabels <- left_join(tmp_box2, df_arc_y, by = c("id" = "grpvar"))
                translabels <- cislabels %>% dplyr::filter(is.na(arc_max_y))
                cislabels <- cislabels %>%
                    dplyr::filter(!is.na(arc_max_y)) %>%
                    distinct(grp, .keep_all = TRUE)
                grid.text(cislabels$txt2, cislabels$arc_center_x + labels.h.offset.cis, cislabels$arc_max_y +
                    annot_size + labels.v.offset.cis,
                rot = rotation.group, gp = font_gpar,
                default.units = "native", just = labels.align
                )
            } else {
                translabels <- tmp_box2
            }
            # TRANS LABELS
            if (nrow(translabels) != 0) {
                translabels$joined_lab <- str_c(translabels$txt1, "_", translabels$txt2)
                grid.text(translabels$joined_lab,
                    translabels$x + labels.h.offset.trans,
                    translabels$yOffset.A + labels.v.offset.trans,
                    rot = rotation.group, gp = font_gpar,
                    default.units = "native", just = labels.align
                )
            }
        } else {
            # plot per-box labels
            grid.text(tmp.A$txt1, tmp.A$x, yOffset.A + annot_size,
                rot = rotation.group, gp = font_gpar,
                default.units = "native", just = labels.align
            )

            grid.text(tmp.B$txt2, tmp.B$x, yOffset.B + annot_size,
                rot = rotation.group, gp = font_gpar,
                default.units = "native", just = labels.align
            )
        }
    }
    popViewport(1)
    # Brought here from AnnotationTrack # EgorS
    ## Finally we set up the image map
    im <- if (!is.null(box)) {
        coords <- as.matrix(box[, c("x1", "y1", "x2", "y2"), drop = FALSE])
        restCols <- setdiff(colnames(box), c("x1", "x2", "y1", "y2", "cx1", "cx2", "cy1", "cy2", "textX", "textY"))
        tags <- lapply(restCols, function(x) {
            tmp <- as.character(box[, x])
            names(tmp) <- rownames(coords)
            tmp
        })
        names(tags) <- restCols
        tags$title <- identifier(GdObject)
        Gviz:::ImageMap(coords = coords, tags = tags)
    } else {
        NULL
    }
    imageMap(GdObject) <- im
    return(invisible(GdObject))
})

#' The default display parameters for a  `DuplexTrack` object
#'
#' `DuplexTrack` inherits from `[Gviz::Annotaiontrack()]` and its Gviz parents.
#'  Most likely, user doesn't need all dioplay pars for the parents, so only
#'  parameters relevant to the `DuplexTrack` are returned by default.
#' @param class  `DuplexTrack` track object
#'   This function allows user to display the
#'   default display parameters for the \code{DuplexTrack} class.
#' @returns list of the default display parameters.
#' @importFrom Gviz availableDisplayPars
#' @export
#' @examples
#' library(InteractionSet)
#' anchor1 <- GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges(
#'         start = c(100, 600, 1100, 1600, 2100),
#'         end = c(200, 700, 1200, 1700, 2200)
#'     ),
#'     strand = "+"
#' )
#' anchor2 <- GRanges(
#'     seqnames = "chr1",
#'     ranges = IRanges(
#'         start = c(300, 800, 1300, 1800, 2300),
#'         end = c(400, 900, 1400, 1900, 2400)
#'     ),
#'     strand = "+"
#' )
#'
#' interactions <- GInteractions(anchor1, anchor2, mode = "strict")
#' gr_region <- range(anchor1, anchor2)
#' a <- DuplexTrack(interactions, gr_region = gr_region, stacking = "dense")
#' availableDisplayPars("DuplexTrack")
#' DuplexDiscoverer::availableDisplayPars(a)
availableDisplayPars <- function(class) {
    if (!is.character(class)) {
        class <- class(class)
    }
    # trick from Ginteraction Interaction to not mask Gviz generic class
    if (class == "DuplexTrack") {
        displaypars <- as(getClassDef(class)@prototype@dp, "list")
        return(displaypars)
    } else {
        return(Gviz::availableDisplayPars(class))
    }
}


#' Plots distributed boxes
#'
#' Non-exported from Gviz, contains logic for calcualting boxes alignment on plot
#'
#' @keywords internal
#' @returns pushes boxes to the viewport
.gv_plotboxes <- function(box, lwd, lty, alpha) {
    if ("transcript" %in% colnames(box)) {
        box <- .handleComposite(box, "box")
        if (nrow(box$box)) {
            grid.rect(box$box$cx2, box$box$cy1,
                width = box$box$cx2 -
                    box$box$cx1, height = box$box$cy2 - box$box$cy1,
                gp = gpar(
                    col = as.character(box$box$col), fill = as.character(box$box$fill),
                    lwd = lwd, lty = lty, alpha = alpha
                ), default.units = "native",
                just = c("right", "bottom")
            )
        }
        if (nrow(box$pols)) {
            grid.polygon(
                x = box$pols$x, y = box$pols$y, id = factor(box$pols$id),
                gp = gpar(
                    col = as.character(box$polpars$col),
                    fill = as.character(box$polpars$fill), lwd = lwd,
                    lty = lty, alpha = alpha
                ), default.units = "native"
            )
        }
    } else {
        grid.rect(box$cx2, box$cy1,
            width = box$cx2 - box$cx1,
            height = box$cy2 - box$cy1, gp = gpar(
                col = as.character(box$col),
                fill = as.character(box$fill), lwd = lwd, lty = lty,
                alpha = alpha
            ), default.units = "native", just = c(
                "right",
                "bottom"
            )
        )
    }
}

#' Takes  of inherited Gviz graphical parameters
#' @importFrom Gviz getScheme
#' @keywords internal
#' @returns set default values inside class upon calling withn class constructor
.gv_updatepars <- function(x, class) {
    current <- getPar(x, hideInternal = FALSE)
    defaults <- getPar(getClass(class)@prototype@dp, hideInternal = FALSE)
    sid <- getOption("Gviz.scheme")
    if (is.null(sid)) {
        scheme <- list()
    } else {
        # use scheme compatibility
        scheme <- Gviz::getScheme(sid)
    }
    # set thisto calculate labels - we use groups as interactions pars
    scheme[["AnnotationTrack"]]["groupAnnotation"] <- "group"
    schemePars <- if (is.null(scheme) || is.null(scheme[[class]])) {
        list()
    } else {
        scheme[[class]]
    }
    if (!is.list(schemePars)) {
        stop(sprintf(
            "Corrupted parameter definition for class '%s' in scheme '%s'",
            class, sid
        ))
    }
    defaults[names(schemePars)] <- schemePars
    missing <- setdiff(names(defaults), names(current))
    if (is.null(getPar(x, ".__appliedScheme")) && length(schemePars)) {
        defaults[[".__appliedScheme"]] <- if (is.null(sid)) {
            NA
        } else {
            sid
        }
        defaults[names(schemePars)] <- schemePars
        missing <- c(union(missing, names(schemePars)), ".__appliedScheme")
    }
    x <- setPar(x, defaults[missing], interactive = FALSE)
    return(x)
}
