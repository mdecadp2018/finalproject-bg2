<?xml version="1.0" encoding="utf-8"?>
<!-- Created by Leo: http://leoeditor.com/leo_toc.html -->
<leo_file xmlns:leo="http://leoeditor.com/namespaces/leo-python-editor/1.1" >
<leo_header file_format="2" tnodes="0" max_tnode_index="0" clone_windows="0"/>
<globals body_outline_ratio="0.5" body_secondary_ratio="0.5">
	<global_window_position top="50" left="50" height="500" width="700"/>
	<global_log_window_position top="0" left="0" height="0" width="0"/>
</globals>
<preferences/>
<find_panel_settings/>
<vnodes>
<v t="amd.20160714225308.1"><vh>@settings</vh>
<v t="amd.20160714225308.2"><vh>@data qt-gui-plugin-style-sheet</vh></v>
<v t="amd.20160714225308.3"><vh>@string initial_split_orientation = horizontal</vh></v>
</v>
<v t="leo.20180702145534.1"><vh>wsgi start/stop</vh>
<v t="leo.20180702142719.1"><vh>@button start wsgi</vh></v>
<v t="leo.20180702142849.1"><vh>@button stop wsgi</vh></v>
</v>
<v t="leo.20180702145555.1"><vh>http start/stop</vh>
<v t="leo.20180702145442.1"><vh>@button start http.server</vh></v>
<v t="leo.20180702145458.1"><vh>@button stop http.server</vh></v>
</v>
<v t="leo.20180703102027.1"><vh>@edit README.md</vh></v>
<v t="amd.20160714225245.2" a="E"><vh>about cmsimfly</vh>
<v t="amd.20160714230137.1"><vh>@edit .gitignore</vh></v>
<v t="amd.20160714230155.1"><vh>@edit LICENSE</vh></v>
<v t="amd.20160715000709.1" a="E"><vh>註解</vh>
<v t="amd.20160715211540.1"><vh>處理 url 與 edit</vh></v>
<v t="amd.20160715213810.1"><vh>必須處理 downloads 與 images 目錄</vh></v>
<v t="leo.20180625155103.1" a="E"><vh>改為 pybean based</vh>
<v t="leo.20180625155310.1"><vh>資料庫設計</vh></v>
</v>
<v t="leo.20180630135634.1" a="E"><vh>新增 generate_pages 功能</vh>
<v t="leo.20180630211407.1"><vh>必須修改頁面內容</vh></v>
<v t="leo.20180712144340.1"><vh>改回手動轉靜態網頁</vh></v>
</v>
</v>
<v t="leo.20180701102916.1" a="E"><vh>2018Summer 改版</vh>
<v t="leo.20180702104314.1"><vh>整合設定檔案</vh></v>
<v t="leo.20180701205822.1"><vh>editorhead 修改</vh></v>
<v t="leo.20180701212101.1"><vh>有關 Cango library</vh></v>
<v t="leo.20180726212905.1"><vh>存檔時移除所有 tag of title</vh></v>
<v t="leo.20180728154823.1"><vh>利用 flake8 查驗格式</vh></v>
</v>
<v t="amd.20160718002238.1"><vh>參考資料</vh></v>
<v t="amd.20160717220157.1"><vh>Flask 開講</vh></v>
</v>
<v t="leo.20180701102750.1"><vh>@edit index.html</vh></v>
<v t="amd.20160714225914.1"><vh>@edit setup.py</vh></v>
<v t="amd.20160714230021.1"><vh>@edit init.py</vh></v>
<v t="amd.20160714230050.1"><vh>@edit nocache.py</vh></v>
<v t="leo.20180625155248.1"><vh>@edit pybean.py</vh></v>
<v t="leo.20180625155728.1"><vh>@edit schema.sql</vh></v>
<v t="amd.20160718130611.1"><vh>@edit .gitignore</vh></v>
<v t="amd.20160717222221.1"><vh>@path static</vh>
<v t="leo.20180718091411.1"><vh>@edit cmsimply.css</vh></v>
<v t="amd.20160717222242.1"><vh>@edit axuploader.js</vh></v>
<v t="amd.20160718125932.1"><vh>@edit viewstl.html</vh></v>
<v t="amd.20160718125943.1"><vh>@edit stlviewer.js</vh></v>
<v t="amd.20160718212341.1"><vh>@edit jsc3d.js</vh></v>
</v>
<v t="amd.20160714230058.1" a="E"><vh>@clean wsgi.py</vh>
<v t="amd.20160714230122.1"><vh>wsgi declarations</vh></v>
</v>
<v t="amd.20160714230833.1" a="E"><vh>編輯 flaskapp.py</vh>
<v t="leo.20181015203721.1" a="E"><vh>@auto flaskapp.py</vh></v>
</v>
</vnodes>
<tnodes>
<t tx="amd.20160714225245.2">we are going to use flask to write an simple cms system based upon cmsimply.

</t>
<t tx="amd.20160714225308.1"></t>
<t tx="amd.20160714225308.2">QSplitter::handle {
    background-color: #CAE1FF; /* lightSteelBlue1 */
}

QStackedWidget {
    /* background-color:lightpink; */
    border-color: red;
    padding: 0px;
    /* border-width: 0px; */
    /* background-color: yellow; */
}

QSplitter {
    border-color: white;
    background-color: white;
    border-style: solid;
}

QTreeWidget {
    /* These apply to the selected item, but not to editing items.*/
    background-color: #ffffec; /* Leo's traditional tree color */
    selection-color: black; /* was white */
    selection-background-color: lightgrey;
    /* font-family: SansSerif; */
    /*font-family: DejaVu Sans Mono;*/
    font-family:YaHei Mono;
    /* 標題字型大小設定 */
    font-size: 22px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal, italic,oblique */
 }

/* Headline edit widgets */
QTreeWidget QLineEdit {
    background-color: cornsilk;
    selection-color: white;
    selection-background-color: blue;
    /*font-family: DejaVu Sans Mono;*/    
    font-family:YaHei Mono;
    /* 沒有特別對應字型大小 */
    font-size: 22px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal, italic,oblique */
}

/* The log panes */
QTextEdit {
    background-color: #f2fdff;
    selection-color: red;
    selection-background-color: blue;
    /* font-family: Courier New; */
    font-family:YaHei Mono;
    /* log font 大小 */
    font-size: 22px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal, italic,oblique */
}

/* The body pane */
QTextEdit#richTextEdit {
    background-color: #fdf5f5; /* A kind of pink. */
    selection-color: white;
    selection-background-color: red;
    /*font-family: DejaVu Sans Mono;*/
    /* font-family: Courier New; */
    font-family:YaHei Mono;
    /* 內文字型大小 */
    font-size: 22px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal,italic,oblique */
}

QLabel {
    font-family:YaHei Mono;
    /* 下方的 Minibuffer 標題字型大小 */
    font-size: 22px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal,italic,oblique */
}

/* Editor labels */
QLineEdit#editorLabel {
    background-color: #ffffec;
    font-family:YaHei Mono;
    /* 沒有直接對應字型大小 */
    font-size: 22px;
    font-weight: normal; /* normal,bold,100,..,900 */
    font-style: normal; /* normal,italic,oblique */
    border: 2px;
    margin: 2px;
}</t>
<t tx="amd.20160714225308.3">horizontal: body pane to the left
vertical: body pane on the botton</t>
<t tx="amd.20160714230058.1">@language python
@tabwidth -4
@others
</t>
<t tx="amd.20160714230122.1">#!/usr/bin/python
# 導入同目錄下的 flaskapp.py
import flaskapp
import ssl

# 即使在近端仍希望以 https 模式下執行
context = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
context.load_cert_chain('localhost.crt', 'localhost.key')

# 取 flaskapp.py 中的 uwsgi 變數設定
uwsgi = flaskapp.uwsgi
if uwsgi:
    # 表示程式在雲端執行
    application = flaskapp.app
else:
    # 表示在近端執行, 以 python3 wsgi.py 執行
    flaskapp.app.run(
        host='127.0.0.1', port=8443, debug=True,
        ssl_context=context
    )
</t>
<t tx="amd.20160714230833.1"></t>
<t tx="amd.20160715000709.1">return redirect(url_for('edit_page')) 這裡的 edit_page 指函式方法

但是若 return redirect("/edit_page') 則必須指向 url

初步完成, 接下來讓 Flask 版本可以正常運作...</t>
<t tx="amd.20160715211540.1">@app.route('/sitemap', defaults={'edit':1})
@app.route('/sitemap/&lt;path:edit&gt;')
def sitemap(edit):
    head, level, page = parse_content()
    directory = render_menu(head, level, page)
    sitemap = render_menu(head, level, page, sitemap=1)
    return set_css()+"&lt;div class='container'&gt;&lt;nav&gt;"+ \
    directory+"&lt;/nav&gt;&lt;section&gt;&lt;h1&gt;Site Map&lt;/h1&gt;"+sitemap+"&lt;/section&gt;&lt;/div&gt;&lt;/body&gt;&lt;/html&gt;"
    
因為 edit=1 針對各頁面標題, 可以從 view mode 轉為 edit mode, 但是使用者也可以在選擇其他功能後, 按下 edit, 這時就必須設法讓各函式可以正常運作

主要關鍵在無 edit 時, 採用內定值, 而隨後無論有幾層 url 變數, 都包含在 path 變數中</t>
<t tx="amd.20160715213810.1"></t>
<t tx="amd.20160717220157.1">get 方法

@app.route('/file_selector', methods=['GET'])
def file_selector():
    if not isAdmin():
        return redirect("/login")
    else:
        type = request.args.get('type')
        page = request.args.get('page', 1)
        item_per_page = request.args.get('item_per_page', 10)
        keyword = request.args.get('keyword')
        if type == "file":
            type = 'downloads'
            return file_lister(download_dir, type, page, item_per_page)
        elif type == "image":
            return file_lister(image_dir, type, page, item_per_page)
            
post 方法

@app.route('/saveConfig', methods=['POST'])
def saveConfig():
    if not isAdmin():
        return redirect("/login")
    site_title = request.form['site_title']
    password = request.form['password']
    password2 = request.form['password2']
    if site_title == None or password == None:
        return error_log("no content to save!")
    old_site_title, old_password = parse_config()
    head, level, page = parse_content()
    directory = render_menu(head, level, page)
    if site_title == None or password == None or password2 != old_password or password == '':
        return set_css()+"&lt;div class='container'&gt;&lt;nav&gt;"+ \
        directory+"&lt;/nav&gt;&lt;section&gt;&lt;h1&gt;Error!&lt;/h1&gt;&lt;a href='/'&gt;Home&lt;/a&gt;&lt;/body&gt;&lt;/html&gt;"
    else:
        if password == password2 and password == old_password:
            hashed_password = old_password
        else:
            hashed_password = hashlib.sha512(password.encode('utf-8')).hexdigest()
        file = open(config_dir+"config", "w", encoding="utf-8")
        file.write("siteTitle:"+site_title+"\npassword:"+hashed_password)
        file.close()
        return set_css()+"&lt;div class='container'&gt;&lt;nav&gt;"+ \
        directory+"&lt;/nav&gt;&lt;section&gt;&lt;h1&gt;config file saved&lt;/h1&gt;&lt;a href='/'&gt;Home&lt;/a&gt;&lt;/body&gt;&lt;/html&gt;"

flask 特殊的 URL 傳值

@app.route('/get_page')
@app.route('/get_page/&lt;heading&gt;', defaults={'edit': 0})
@app.route('/get_page/&lt;heading&gt;/&lt;int:edit&gt;')
def get_page(heading, edit):
    head, level, page = parse_content()
    directory = render_menu(head, level, page)
    if heading == None:
        heading = head[0]
    # 因為同一 heading 可能有多頁, 因此不可使用 head.index(heading) 搜尋 page_order
    page_order_list, page_content_list = search_content(head, page, heading)
    return_content = ""
    pagedata = ""
    outstring = ""
    outstring_duplicate = ""
    pagedata_duplicate = ""
    outstring_list = []
    for i in range(len(page_order_list)):
        #page_order = head.index(heading)
        page_order = page_order_list[i]
        if page_order == 0:
            last_page = ""
        else:
            last_page = head[page_order-1]+" &lt;&lt; &lt;a href='/get_page/"+head[page_order-1]+"'&gt;Previous&lt;/a&gt;"
        if page_order == len(head) - 1:
            # no next page
            next_page = ""
        else:
            next_page = "&lt;a href='/get_page/"+head[page_order+1]+"'&gt;Next&lt;/a&gt; &gt;&gt; "+ head[page_order+1]
        if len(page_order_list) &gt; 1:
            return_content += last_page+" "+next_page+"&lt;br /&gt;&lt;h1&gt;"+heading+"&lt;/h1&gt;"+page_content_list[i]+"&lt;br /&gt;"+last_page+" "+next_page+"&lt;br /&gt;&lt;hr&gt;"
            pagedata_duplicate = "&lt;h"+level[page_order]+"&gt;"+heading+"&lt;/h"+level[page_order]+"&gt;"+page_content_list[i]
            outstring_list.append(last_page+" "+next_page+"&lt;br /&gt;"+ tinymce_editor(directory, cgi.escape(pagedata_duplicate), page_order))
        else:
            return_content += last_page+" "+next_page+"&lt;br /&gt;&lt;h1&gt;"+heading+"&lt;/h1&gt;"+page_content_list[i]+"&lt;br /&gt;"+last_page+" "+next_page
            
        pagedata += "&lt;h"+level[page_order]+"&gt;"+heading+"&lt;/h"+level[page_order]+"&gt;"+page_content_list[i]
        outstring += last_page+" "+next_page+"&lt;br /&gt;"+ tinymce_editor(directory, cgi.escape(pagedata), page_order)
    
    # edit=0 for viewpage
    if edit == 0:
        return set_css()+"&lt;div class='container'&gt;&lt;nav&gt;"+ \
        directory+"&lt;/nav&gt;&lt;section&gt;"+return_content+"&lt;/section&gt;&lt;/div&gt;&lt;/body&gt;&lt;/html&gt;"
    # enter edit mode
    else:
        # check if administrator
        if not isAdmin():
            redirect(url_for('login'))
        else:
            if len(page_order_list) &gt; 1:
                # 若碰到重複頁面頁印, 且要求編輯, 則導向 edit_page
                #return redirect("/edit_page")
                for i in range(len(page_order_list)):
                    outstring_duplicate += outstring_list[i]+"&lt;br /&gt;&lt;hr&gt;"
                return outstring_duplicate
            else:
            #pagedata = "&lt;h"+level[page_order]+"&gt;"+heading+"&lt;/h"+level[page_order]+"&gt;"+search_content(head, page, heading)
            #outstring = last_page+" "+next_page+"&lt;br /&gt;"+ tinymce_editor(directory, cgi.escape(pagedata), page_order)
                return outstring</t>
<t tx="amd.20160717222221.1"></t>
<t tx="amd.20160718002238.1">有關內容解碼

http://stackoverflow.com/questions/12247681/tinymce-is-adding-nbsp-instead-of-the-space-when-using-the-word-paste

現在標題中有空白時會使用 + 但是正確的編碼應該是 %20</t>
<t tx="leo.20180625155103.1">2018.06.25 啟動

使用 SQLite 資料庫儲存檔案, 而且可以轉為純 html 文件

2018.07.25 暫時關掉 pybean

</t>
<t tx="leo.20180625155310.1">CMSimfly 所需欄位

id - 好像是 pybean 對應資料庫必須有的 index 欄位
follow - default 為 0 表示為 level1, 其餘按照 follow 順利遞迴列出, 理論上可以無限階
title - 文章標題
content - 文章內容
memo - 額外註解欄位

</t>
<t tx="leo.20180630135634.1">可以將 CMSimfly 轉為 html pages 
其中以 index.html 為起始點
結合 downloads, images, static 與 content 目錄中的子頁面, 可以放入 Github Pages 或 Fossil SCM doc 目錄中, 成為靜態網站內容.

</t>
<t tx="leo.20180630211407.1">圖片的連結必須改為靜態模式, 而非目前的動態模式

整體依附在 OpenShift 的部分必須刪除?

可能在 generate_pages  時針對 downloads 與 images 目錄中的引用,  將動態格式轉為靜態格式</t>
<t tx="leo.20180701102916.1">去除 OpenShift 相關程式碼

增加轉為靜態網頁格式功能

整體目錄架構修改

static 目錄存放各種 css 設定與 Javascript 檔案

downloads 目錄存放使用者在動態模式下上傳的文件檔

images 目錄存放使用者在動態模式下上傳的影像檔

content 目錄存放靜態 html 檔案

config 目錄存放 config 設定檔案與 content.htm - 使用者若在區域網路中執行, 可以將 config 檔案設定到其他目錄中存放, 目錄是可不將 config 檔案新增提交到版本系統

wsgi.py 與其他 Flask based 程式碼

index.html 在純 WWW 模式下 redirect 到 content/index.html

靜態網頁引用 static, downloads, images 目錄內容時, 以相對目錄指定 (./../)

CMSimfly 之後將可以在動態或靜態模式下執行運作</t>
<t tx="leo.20180701205822.1">原來設定

def editorhead():
    return '''
    &lt;br /&gt;
&lt;script src="//cdn.tinymce.com/4/tinymce.min.js"&gt;&lt;/script&gt;
&lt;script src="/static/tinymce4/tinymce/plugins/sh4tinymce/plugin.min.js"&gt;&lt;/script&gt;
&lt;link rel = "stylesheet" href = "/static/tinymce4/tinymce/plugins/sh4tinymce/style/style.css"&gt;
&lt;script&gt;
tinymce.init({
  selector: "textarea",
  height: 500,
  element_format : "xhtml",
  language : "en",
  plugins: [
    'advlist autolink lists link image charmap print preview hr anchor pagebreak',
    'searchreplace wordcount visualblocks visualchars code fullscreen',
    'insertdatetime media nonbreaking save table contextmenu directionality',
    'emoticons template paste textcolor colorpicker textpattern imagetools sh4tinymce'
  ],
  toolbar1: 'insertfile save undo redo | styleselect | bold italic | alignleft aligncenter alignright alignjustify | bullist numlist outdent indent',
  toolbar2: 'link image | print preview media | forecolor backcolor emoticons | code sh4tinymce',
  relative_urls: false,
  toolbar_items_size: 'small',
  file_picker_callback: function(callback, value, meta) {
        cmsFilePicker(callback, value, meta);
    },
  templates: [
    { title: 'Test template 1', content: 'Test 1' },
    { title: 'Test template 2', content: 'Test 2' }
  ],
  content_css: [
    '//fonts.googleapis.com/css?family=Lato:300,300i,400,400i',
    '//www.tinymce.com/css/codepen.min.css'
  ]
});

function cmsFilePicker(callback, value, meta) {
    tinymce.activeEditor.windowManager.open({
        title: 'Uploaded File Browser',
        url: '/file_selector?type=' + meta.filetype,
        width: 800,
        height: 550,
    }, {
        oninsert: function (url, objVals) {
            callback(url, objVals);
        }
    });
};
&lt;/script&gt;
'''


修改後

def editorhead():
    return '''
    &lt;br /&gt;
&lt;script src="//cdn.tinymce.com/4/tinymce.min.js"&gt;&lt;/script&gt;
&lt;script src="/static/tinymce4/tinymce/plugins/sh4tinymce/plugin.min.js"&gt;&lt;/script&gt;
&lt;link rel = "stylesheet" href = "/static/tinymce4/tinymce/plugins/sh4tinymce/style/style.css"&gt;
&lt;script&gt;
tinymce.init({
  selector: "textarea",
  height: 500,
  element_format : "html",
  language : "en",
  valid_elements : '*[*]',
  extended_valid_elements: "script[language|type|src]",
  plugins: [
    'advlist autolink lists link image charmap print preview hr anchor pagebreak',
    'searchreplace wordcount visualblocks visualchars code fullscreen',
    'insertdatetime media nonbreaking save table contextmenu directionality',
    'emoticons template paste textcolor colorpicker textpattern imagetools sh4tinymce'
  ],
  toolbar1: 'insertfile save undo redo | styleselect | bold italic | alignleft aligncenter alignright alignjustify | bullist numlist outdent indent',
  toolbar2: 'link image | print preview media | forecolor backcolor emoticons | code sh4tinymce',
  relative_urls: false,
  toolbar_items_size: 'small',
  file_picker_callback: function(callback, value, meta) {
        cmsFilePicker(callback, value, meta);
    },
  templates: [
    { title: 'Test template 1', content: 'Test 1' },
    { title: 'Test template 2', content: 'Test 2' }
  ],
  content_css: [
    '//fonts.googleapis.com/css?family=Lato:300,300i,400,400i',
    '//www.tinymce.com/css/codepen.min.css'
  ]
});

function cmsFilePicker(callback, value, meta) {
    tinymce.activeEditor.windowManager.open({
        title: 'Uploaded File Browser',
        url: '/file_selector?type=' + meta.filetype,
        width: 800,
        height: 550,
    }, {
        oninsert: function (url, objVals) {
            callback(url, objVals);
        }
    });
};
&lt;/script&gt;
'''


</t>
<t tx="leo.20180701212101.1">&lt;script type="text/javascript" src="/static/Cango-13v08-min.js"&gt;&lt;/script&gt;
&lt;script type="text/javascript" src="/static/CangoAxes-4v01-min.js"&gt;&lt;/script&gt;
&lt;script type="text/javascript" src="/static/gearUtils-05.js"&gt;&lt;/script&gt;

必須直接放在 syntaxhighlight 與 syntaxhighlight2 中, 因為若直接放入各頁面, 會導致只有動態系統可以正常運作, 而靜態頁面無法擷取, 因為動態採絕對目錄取 static 內容, 而靜態頁面則必須採相對目錄擷取.

</t>
<t tx="leo.20180702104314.1">目前以 uwsgi 設定判斷是在一般模式執行 (python wsgi.py) 或 uwsgi 模式執行 (http://lab.kmol.info/blog/cmsimfly-and-fossil-on-same-machine.html)

uwsgi = False

則在一般模式下執行

uwsgi = True

則在 uwsgi 模式下執行

原始的 uwsgi 設定在 init.py 中的 Init 類別中, 設為 class variable, 也就是所謂的 static variable

然後 flaskapp.py 導入 init.py 後, 將 uwsgi 變數設為 init.Init.uwsgi

最後 wsgi.py 又導入 flaskapp.py, 再將其內的 uwsgi 設為 flaskapp.uwsgi, 就可以輾轉取出  init.py 中的設定</t>
<t tx="leo.20180702142719.1">@language python
bpm = g.app.backgroundProcessManager
bpm.start_process(c, command='python wsgi.py', kind='wsgi', shell=True)


</t>
<t tx="leo.20180702142849.1">@language python
bpm = g.app.backgroundProcessManager
bpm.kill(kind="wsgi") 

</t>
<t tx="leo.20180702145442.1">@language python
bpm = g.app.backgroundProcessManager
bpm.start_process(c, command='python -m http.server', kind='http', shell=True)


</t>
<t tx="leo.20180702145458.1">@language python
bpm = g.app.backgroundProcessManager
bpm.kill(kind="http") 

</t>
<t tx="leo.20180702145534.1"></t>
<t tx="leo.20180702145555.1"></t>
<t tx="leo.20180712144340.1">因為每次 generate_pages 必須將 content 目錄中的所有檔案刪除後, 再從 content.htm 轉為靜態網頁, 因此從負擔較重的每次存檔都自動轉檔, 改回必要時再以手動呼叫 generate_pages 進行轉檔.</t>
<t tx="leo.20180726212905.1">存檔時移除所有 tag of title

這是一個 bug, 目前若 title 有 tag, 則無法透過 get_page 辨識, 因此必須在存檔前移除所有 tags

bug 180726</t>
<t tx="leo.20180728154823.1">將程式碼逐一修改為 flake8 查驗無誤</t>
</tnodes>
</leo_file>
