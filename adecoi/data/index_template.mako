<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="https://raw.githubusercontent.com/cov-ert/civet/master/docs/virus.svg">

    <title>ADE3 ${year} Sequencing Results</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5zCYotlSAcp1+c8xmyTe9GYg1l9a69psu" crossorigin="anonymous">
    <script src="https://code.jquery.com/jquery-1.12.4.min.js" integrity="sha384-nvAa0+6Qg9clwYCGGPpDQLVpLNn0fRaROjHqs13t4Ggj3Ez50XnGQqc/r8MhnRDZ" crossorigin="anonymous"></script>
    <link href="https://cdn.datatables.net/1.10.25/css/jquery.dataTables.min.css" rel="stylesheet" type="text/css"/>
    <script src="https://cdn.datatables.net/1.10.25/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/dataTables.buttons.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.html5.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.print.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha384-aJ21OjlMXNL5UyIl/XNwTMqvzeRMZH2w8c5cRVpzpU8Y5bApTppSuUkhZXN0VxHd" crossorigin="anonymous"></script>
    <script src="https://cdn.jsdelivr.net/gh/rambaut/figtree.js@9880/dist/figtree.umd.js"></script>
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <script src="https://sharonchoong.github.io/svg-exportJS/svg-export.min.js"></script>
    <script src="https://unpkg.com/canvg/lib/umd.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/pdfkit/js/pdfkit.min.js"></script>
    <script src="https://github.com/devongovett/blob-stream/releases/download/v0.1.3/blob-stream.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/svg-to-pdfkit/source.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega@5.16.0"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@4.15.0"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6.11.1"></script>

    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <style>
    body {
        padding-top: 50px;
        font-family: "ArialNova-Light", "HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
    }

    table text {
        font-family: "ArialNova-Light", "HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
    }

    header {
        display: block;
        text-align: right;
    }

    svg {
        width: 90%;
        height: auto;
    }

    .row {
        display: flex;
    }

    .column {
        padding: 10px;
        flex: 50%;
    }

    .accordion {
        background-color: #eee;
        color: #444;
        cursor: pointer;
        padding: 13px;
        width: 100%;
        border: none;
        text-align: left;
        outline: none;
        transition: 0.4s;
    }

    .active, .accordion:hover {
        background-color: #419595;
        color: white;
    }

    .accordion:after {
        content: '\002B';
        color: white;
        font-weight: bold;
        float: right;
        margin-left: 5px;
    }

    .active:after {
        content: "\2212";
    }

    .panel {
        padding: 0 13px;
        background-color: white;
        max-height: 0;
        overflow: hidden;
        transition: max-height 0.2s ease-out;
    }

    .center {
        display: block;
        margin-left: auto;
        margin-right: auto;
        width: 50%;
    }

    .node-background {
        fill: dimgrey;
        stroke: dimgrey;
    }

    .node circle {
        stroke-width: 0;
        cursor: pointer;
        /* fill:#7178bc; */
        stroke: dimgrey;
    }

    .node circle.selected {
        stroke-width: 0;
        cursor: pointer;
        fill: #419595;
        stroke: dimgrey;
    }

    .node-background.query_boolean-True {
        stroke: #419595;
    }

    .node.query_boolean-True circle {
        stroke: #419595;
    }

    .node.query_boolean-True circle.selected {
        stroke: #419595;
    }

    .node-background.query_boolean-True circle.selected {
        stroke: #419595;
    }

    .node.query_boolean-True.hovered circle {
        stroke: #419595;
    }

    .node rect {
        stroke-width: 2;
        fill: #419595;
        stroke: dimgrey;
    }

    .svg-tooltip {
        background: rgba(69, 77, 93, .9);
        border-radius: .1rem;
        color: #fff;
        display: block;
        font-size: 14px;
        max-width: 320px;
        padding: .2rem .4rem;
        position: absolute;
        text-overflow: ellipsis;
        white-space: pre;
        z-index: 300;
        visibility: hidden;
    }

    .tooltip-header {
        font-size: 1.3em;
    }

    .tooltip-key {
        font-weight: bold;
    }

    .branch path {
        stroke-width: 2;
        stroke: dimgrey;
        stroke-linejoin: round;
        cursor: pointer;
    }

    .branch.hovered path {
        stroke-width: 4;
        stroke: dimgrey;
    }

    .node.hovered circle {
        stroke-width: 5;
        stroke: dimgrey
    }

    .node text {
        font-family: "ArialNova-Light", "HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
        font-weight: 300;
        font-size: 0.9em;
    }/* .starter-template {
            padding: 40px 15px;
            text-align: left;
          } */

    .dataTables_wrapper.no-footer .dataTables_scrollBody {
        border-top: 1px solid rgb(148, 148, 148);
        border-bottom: none;
    }

    .svg-icon {
        display: inline-flex;
        align-self: center;
    }

    h3 {
        font-size: 1em;
    }

    #toTopBtn {
        position: fixed;
        bottom: 26px;
        right: 39px;
        z-index: 98;
        padding: 21px;
        background-color: #419595
    }

    .js .cd-top--fade-out {
        opacity: .5
    }

    .js .cd-top--is-visible {
        visibility: visible;
        opacity: 1
    }

    .js .cd-top {
        visibility: hidden;
        opacity: 0;
        transition: opacity .3s, visibility .3s, background-color .3s
    }

    .cd-top {
        position: fixed;
        bottom: 20px;
        bottom: var(--cd-back-to-top-margin);
        right: 20px;
        right: var(--cd-back-to-top-margin);
        display: inline-block;
        height: 40px;
        height: var(--cd-back-to-top-size);
        width: 40px;
        width: var(--cd-back-to-top-size);
        box-shadow: 0 0 10px rgba(0, 0, 0, .05) !important;
        background: url(https://res.cloudinary.com/dxfq3iotg/image/upload/v1571057658/cd-top-arrow.svg) no-repeat center 50%;
        background-color: #419595;
        background-color: hsla(var(--cd-color-3-h), var(--cd-color-3-s), var(--cd-color-3-l), 0.8)
    }

    .slidecontainer {
        width: 100%;
    }

    .colourSelect {
        background: #eee;
        border-radius: 5px;
        padding: 4px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
    }

    .slider {
        -webkit-appearance: none;
        width: 100%;
        height: 15px;
        background: #d3d3d3;
        border-radius: 5px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
        -webkit-transition: .2s;
        transition: opacity .2s;
    }

    .slider:hover {
        opacity: 1;
    }

    .slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 25px;
        height: 25px;
        border-radius: 50%;
        background: #419595;
        stroke: dimgrey;
        cursor: pointer;
    }

    .slider::-moz-range-thumb {
        width: 25px;
        height: 25px;
        border-radius: 50%;
        stroke: dimgrey;
        background: #419595;
        cursor: pointer;
    }

    .tree-container {
        max-height: 1000px;
        overflow: scroll;
    }

    .label {
        display: none;
    }

    .label.show {
        display: inline;
    }

    .node.hovered .label {
        display: inline;
    }

    div.sticky {
        position: -webkit-sticky;
        /* Safari */
        position: sticky;
        top: 0;
    }

    .searchbar {
        border-style: solid;
        border-color: lightgrey;
        border-radius: 5px;
        float: right
    }

    @media print {
        .tree-container {
            max-height: none;
            overflow: visible;
        }

        .slider-block {
            display: none;
        }

        .container {
            padding-right: 1.5cm;
            padding-left: 1.5cm;
            padding-bottom: 1.5cm;
            margin: 1cm;
            min-width: 2200px;
            font-size: 2.5vw;
        }

        .searchbar {
            display: none;
        }

        h3 {
            font-size: 2.5vw;
        }

        h2 {
            font-size: 4vw;
            padding: 1cm;
        }

        h1 {
            font-size: 5vw;
        }

        .command-block {
            display: none;
        }

        pre {
            display: none;
        }

        .civet-logo {
            width: 2cm;
            height: 2cm;
        }

        .tree_svg {
            width: 1200px
        }

        .page-footer {
            display: none;
        }

        .civet-header {
            text-align: left;
        }

        .content-block, p {
            page-break-inside: avoid;
        }
    }

    @media screen and (prefers-color-scheme: dark) {
        body {
            background-color: #17141F;
            color: #F2E7DC;
            opacity: 0.95;
        }

        img {
            filter: brightness(.8) contrast(1.2);
        }

        .component {
            background-color: #2C2640;
        }

        .table-striped > tbody > tr:nth-child(odd) {
            background-color: #2C2640;
            /* border-top-color: #3B325B; */
            /* border-color: #2C2640; */
            opacity: 0.9;
        }

        .table {
            border-top: 0px;
            /* border-color: #3B325B; */
            background-color: #17141F;
        }

        .table > tbody > tr > td, .table > tbody > tr > th, .table > tfoot > tr > td, .table > tfoot > tr > th, .table > thead > tr > td, .table > thead > tr > th {
            border-top: none;
        }

        .accordion {
            background-color: #2C2640;
            color: #F2E7DC;
            cursor: pointer;
            padding: 12px;
            width: 100%;
            border: none;
            text-align: left;
            outline: none;
            transition: 0.4s;
        }

        .active, .accordion:hover {
            background-color: #17141F;
        }

        .accordion:after {
            content: '\002B';
            color: #F2E7DC;
            font-weight: bold;
            float: right;
            margin-left: 5px;
        }

        .active:after {
            content: "\2212";
        }

        .panel {
            padding: 0 12px;
            background-color: #2C2640;
            max-height: 0;
            overflow: hidden;
            transition: max-height 0.2s ease-out;
        }

        pre {
            background-color: #3B325B;
            color: #F2E7DC;
            border: none;
            opacity: 0.8;
        }

        .searchbar {
            background-color: #3B325B;
            color: #F2E7DC;
            border-style: none;
            opacity: 0.8;
        }

        .slider {
            background: #F2E7DC;
            stroke: #F2E7DC;
        }

        .slider::-webkit-slider-thumb {
            background: #5F9C82;
            fill: #5F9C82;
            stroke: #F2E7DC;
        }

        .slider::-moz-range-thumb {
            stroke: #F2E7DC;
            background: #5F9C82;
            fill: #5F9C82;
        }

        .node-background {
            fill: #F2E7DC;
            stroke: #F2E7DC;
            opacity: 0.85;
        }

        .node circle {
            stroke-width: 0;
            cursor: pointer;
            fill: #5F9C82;
            stroke: #F2E7DC;
        }

        .node circle.selected {
            stroke-width: 0;
            cursor: pointer;
            fill: #E27E7E;
            stroke: #F2E7DC;
            opacity: 1;
        }

        .node rect {
            stroke-width: 2;
            fill: #E27E7E;
            stroke: #F2E7DC;
        }

        .svg-tooltip {
            background: rgba(69, 77, 93, .9);
            color: #F2E7DC;
        }

        .branch path {
            stroke: #F2E7DC;
            opacity: 0.85;
        }

        .branch.hovered path {
            stroke: #F2E7DC;
            opacity: 1;
        }

        .node.hovered circle {
            stroke: #F2E7DC;
            opacity: 1;
        }

        .node text {
            font-family: "ArialNova-Light", "HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
            font-weight: 300;
            font-size: 0.9em;
            color: #F2E7DC;
            fill: #F2E7DC;
        }

        .scale-bar line {
            stroke: #F2E7DC;
        }

        .scale-bar text {
            fill: #F2E7DC;
            color: #F2E7DC;
        }
    }
    </style>

</head>

<body>
    <div class="container">
        <a href="#" id="toTopBtn" class="cd-top text-replace js-cd-top cd-top--is-visible cd-top--fade-out" data-abc="true"></a>
        <div>
            <header class="civet-header">

                            ADE3 | 
                            
                <small class="text-muted">Animal diversity and evolution 3</small>
                <hr>
            </header>

            <h1>
                ADE3 ${year} Sequencing Results
                            
                <small class="text-muted" style="color:#419595">${date}</small>
            </h1>
            <hr>
            <br>
        </div>
        <br>
        <table class="display nowrap" id="myTable">
            <thead>
                <tr>
                    <th style="width:30%;">Name</th>
                    <th style="width:10%;">ID</th>
                    <th style="width:20%;"></th>
                    <th style="width:20%;"></th>
                    <th style="width:20%;"></th>
                </tr>
            </thead>
                <tbody>
                    % for barcode_data in data_for_report:
                        <% name,barcode = barcode_data %>
                        <tr><td>${name}</td>
                            <td>${barcode}</td>
                            <td><a href="https://uofedinburghade3.github.io/ADE_COI/ADE3_${year}/barcode${barcode}/analysis_report.html">sequencing report</a></td>
                            <td><a href="https://uofedinburghade3.github.io/ADE_COI/ADE3_${year}/barcode${barcode}/krona.filtered.html">read classification</a></td>
                            <td><a href="https://github.com/UofEdinburghADE3/ADE_COI/tree/main/ADE3_${year}/barcode${barcode}/consensus_sequences">consensus sequences</a></td></tr>
                    % endfor
                    </tbody>
        </table>

        <script type="text/javascript">
        $(document).ready(function() {
            var table = $('#myTable').DataTable({
                "paging": false,
                "border-bottom": false,
                dom: 'frtip',
                buttons: ["copy", "csv", "print"]
            });
            table.buttons().container().appendTo($('#tableExportID'));
            $('a.toggle-vis').on('click', function(e) {
                e.preventDefault();

                // Get the column API object
                var column = table.column($(this).attr('data-column'));

                // Toggle the visibility
                column.visible(!column.visible());
            });

        });
        </script>

        <footer class="page-footer">
            <div class="container-fluid text-right text-md-right">
                <hr>
                <div class="row"></div>

                <div class="col-sm-11" style="text-align: right;">

                            ADE3 | 
                    <small class="text-muted">Animal diversity and evolution</small>
                    <br>
                    <small class="text-muted">GNU General Public License v3.0</small>
                </div>

                <br>
                <br>
                </p>
            </div>
        </footer>
    </div>
</body>
</html>
