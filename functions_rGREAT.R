library(rGREAT)
library(rtracklayer)
library(xlsx)

get_type = function(str){
  out = "string"
  if(suppressWarnings(sum(is.na(as.numeric(str))) == 0)){
    out = "numeric"
    if(all(as.numeric(str) == as.integer(str))){
      out = "integer"
    }
  }
  return(out)
}

df2sheet = function(df, sheetName, wb){
  mat = as.matrix(df)
  if(nrow(mat) > 100) mat = mat[1:100,]
  sheet = createSheet(wb = wb, sheetName = sheetName)
  
  for(i in 1:ncol(mat)){
    cb = CellBlock(sheet = sheet, startRow = 1, startColumn = 1, noRows = 1, noColumns = ncol(mat))
    CB.setRowData(cb, colnames(mat), rowIndex = 1)
    
    cb = CellBlock(sheet = sheet, startRow = 2, startColumn = i, noRows = nrow(mat), noColumns = 1)
    col_type = get_type(mat[,i])
    vals = switch (col_type,
                   string = mat[,i],
                   numeric = as.numeric(mat[,i]),
                   integer = as.integer(mat[,i])
    )
    style = switch (col_type,
                    string = NULL,
                    numeric = CellStyle(wb) + DataFormat("0.00E+00"),
                    integer = CellStyle(wb) + DataFormat("0")
    )
    if(col_type == "numeric"){
      if(median(vals) > .1 && median(vals) < 10000){
        style = CellStyle(wb) + DataFormat("#,##0.00")
      }
    }
    # print(col_type)
    # print(class(vals))
    CB.setColData(cellBlock = cb, x = vals, colIndex = 1, colStyle = style)
  }
  
  for(i in 1:ncol(mat)){
    autoSizeColumn(sheet = sheet, colIndex = i)
  }
}

wrap_great = function(subject38, query38, subject_name, query_name, with = T){
  index = 1
  if(!with) index = 2
  
  venn_in = list(subject38,  query38)
  names(venn_in) = c(subject_name, query_name)
  venn38 = venn_peaks(venn_in)
  
  
  ch = import.chain("H:/liftOver_chains/hg38ToHg19.over.chain")
  subject19 = unlist(liftOver(subject38, ch))
  query19 = unlist(liftOver(query38, ch))
  
  subject19  = reduce(subject19)
  
  query19  = reduce(query19)
  
  venn_in = list(subject19,  query19)
  names(venn_in) = c(subject_name, query_name)
  venn = venn_peaks(venn_in)
  
  
  fg19 = venn[[index]]
  bg19 = subject19
  
  job = submitGreatJob(gr = fg19, bg = bg19, bgChoice = "data",  species = "hg19", version = "3.0.0", request_interval = 30)
  wb = createWorkbook()
  
  for(cat in rGREAT::availableCategories(job)){
    if(cat == "Gene Expression") next #Gene Expression doesn't work
    tables = getEnrichmentTables(job, category = cat)
    for(tab_name in names(tables)){
      print(tab_name)
      mat = tables[[tab_name]]
      df2sheet(df = mat, sheetName = tab_name, wb = wb)
    }
  }
  
  saveWorkbook(wb, paste0(subject_name, ifelse(with, " with ", " no "), query_name, ".xlsx"))
  
  return(list(job = job, subject = subject38, bound_query = venn38[[1]], absent_query = venn38[[2]] ))
}
