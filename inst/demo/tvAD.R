library(TrenaViz)
PORT <- 5899
tv <- TrenaViz("TrenaProjectAD")
runApp(createApp(tv, port=PORT))
url <- sprintf("http://0.0.0.0:%d", PORT)
later(function(){browseURL(url)}, 2)
