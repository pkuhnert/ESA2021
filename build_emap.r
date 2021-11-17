build_emap <- function (data, pdflist = NULL, geoData = NULL, id = NULL, key_label, 
    palette = "YlOrRd", size = NULL, border = NULL) 
{
    if (is.null(geoData)) {
        l1 <- match("long", names(data))
        l2 <- match("lat", names(data))
        if (any(is.na(c(l1, l2)))) 
            stop("There must be coordinates in data columns named long and lat.\n")
    }
    if (!is.null(geoData) & is.null(id)) {
        stop("Missing id. Must be a common column shared by data and geoData.")
    }
    if (!is.null(border)) {
        m <- (border %in% c("county", "france", "italy", "nz", 
            "state", "usa", "world"))
        if (length(m) == 0) 
            stop("Border name not recognised. Must be one of county, france, italy, nz, state, usa or world \n\n           (see documentation for map_data function in ggplot2 for more information)")
        else bord <- map_data(border)
    }
    else {
        long <- numeric(0)
        lat = numeric(0)
        bord <- data.frame(long = long, lat = lat, group = numeric(0))
    }
    if (!is.null(size)) {
        s <- seq(1, 20, by = 1)
        if (!(size %in% s) & !is.null(size)) 
            stop("Size not recognised. Must be an integer between 1 and 20.")
    }
    if (!(palette %in% c("YlOrBr", "YlOrRd", "YlGnBu", "PuBuGn"))) 
        stop("Palette name not recognised. Must be one of YlOrBr, YlOrRd, YlGnBu or PuBuGn, \n\n         (see documentation for scale_fill_distiller in ggplot2 for more information).\n")
    if (!is.null(geoData)) {
      
        geoData@data <- geoData@data %>% dplyr::mutate_if(is.factor, 
            as.character)
        geoData@data <- left_join(geoData@data, data, by = id)
        geoData@data$id <- rownames(geoData@data)
        region_coord <- sptable(geoData, region = "id")
        region_coord <- plyr::rename(region_coord, c(object_ = "id", 
            x_ = "long", y_ = "lat", branch_ = "group"))
        output_data <- plyr::join(region_coord, geoData@data, 
            by = "id")
        bbox <- make_bbox(lat = lat, lon = long, data = output_data)
    }
    else {
        output_data <- data
        bbox <- make_bbox(lat = lat, lon = long, data = output_data)
    }
    if (is.null(pdflist)) {
        id <- match(names(data)[1:3], names(output_data))
        names(output_data)[id] <- c("estimate", "error", "pr_exc")
    }
    else {
        
        warning("Ensure the pdf you select is suitable for your data. See ??build_emap for examples of good and bad distribution choices.\n")
        id <- match(names(data)[1:2], names(output_data))
        names(output_data)[id] <- c("estimate", "error")
        estimate <- output_data$estimate
        error <- output_data$error
        args_call <- eval(do.call("substitute", list(pdflist$args, 
            list(estimate, error))))
        output_data$pr_exc <- pdist(pname = pdflist$dist, th = pdflist$th, 
            args = args_call)
    }
    p <- list(output_data = output_data, bord = bord, bbox = bbox, 
        key_label = key_label, palette = palette)
    oldClass(p) <- c("emap", class(p))
    p
}
