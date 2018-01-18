class NotFoundError(Exception):
    '''Raised when node was not found'''


class IntegrityError(Exception):
    '''Raised when tree is broken'''

class DBClade(TreeElement, TreeMixin):
    def __init__(self, adaptor, taxon_id, parent_id, ncbi_id, name, left_val, right_val, branch_length=None, name=None, clades=None,
            confidence=None, color=None, width=None):
        """Define parameters for the Clade tree."""
        self.branch_length = branch_length
        self.name = name
        self.clades = clades or []
        self.confidence = confidence
        self.color = color
        self.width = width
        self.ncbi_id = ncbi_id
        self.name = name
        # the following are private used for SQL queries
        self._adaptor = adaptor
        self._id = taxon_id
        self._parent_id = parent_id
        self._left_val = left_val
        self._right_val = right_val

    def is_terminal(self):
        return self._right_val - self._left_val == 1


class TaxonNode(object):
    def __init__(self, adaptor, taxon_id, parent_id, ncbi_id, name=None, rank=None, genetic_code=None, mito_genetic_code=None, left_val=None, right_val=None, **kwargs):
        ''' Initialize a specific node in the tree.

            Normally a user would not need to construct these objects
            by themselves, instead they should be constructed automatically
            using the find() method of TaxonTree
        '''
        self.ncbi_id = ncbi_id
        self.name = name
        self.rank = rank
        self.genetic_code = genetic_code
        self.mito_genetic_code = mito_genetic_code

        self._adaptor = adaptor
        self._id = taxon_id
        self._parent_id = parent_id
        self._left_val = left_val
        self._right_val = right_val
        for key, value in kwargs.items():
            setattr(self, key, value)

    def is_terminal(self):
        return self._right_val - self._left_val == 1


class TaxonTree(object):

    def __init__(self, adaptor):
        self.adaptor = adaptor

    def _update_left_right_taxon_values(self, left_value, width):
        """update the left and right values in the table
        """
        # Due to the UNIQUE constraint on the left and right values in the taxon
        # table we cannot simply update them through an SQL statement as we risk
        # colliding values. Instead we must select all of the rows that we want to
        # update, modify the values in python and then update the rows
        # self.adaptor.execute("UPDATE taxon SET right_value = right_value + 2 "
        #                      "WHERE right_value >= %s", (left_value,))
        # self.adaptor.execute("UPDATE taxon SET left_value = left_value + 2 "
        #                      "WHERE left_value > %s", (left_value,))

        sql = "SELECT left_value, right_value, taxon_id FROM taxon " \
              "WHERE right_value >= %s OR left_value > %s"
        rows = self.adaptor.execute_and_fetchall(sql, (left_value, left_value))

        right_rows = []
        left_rows = []
        for row in rows:
            new_right = row[1]
            new_left = row[0]
            if new_right >= left_value:
                new_right += width

            if new_left > left_value:
                new_left += width

            right_rows.append((new_right, row[2]))
            left_rows.append((new_left, row[2]))

        # sort the rows based on the value from largest to smallest
        # should ensure no overlaps
        if width < 0:
            rev = False
        else:
            rev = True
        right_rows = sorted(right_rows, key=lambda x: x[0], reverse=rev)
        left_rows = sorted(left_rows, key=lambda x: x[0], reverse=rev)
        for row in left_rows:
            self.adaptor.execute("UPDATE taxon SET left_value = %s"
                                 " WHERE taxon_id = %s", row)

        for row in right_rows:
            self.adaptor.execute("UPDATE taxon SET right_value = %s"
                                 " WHERE taxon_id = %s", row)

    def _make_node(self, _id):
        taxon_sql = '''SELECT taxon_id,
                        ncbi_taxon_id,
                        parent_taxon_id,
                        node_rank,
                        genetic_code,
                        mito_genetic_code,
                        left_value,
                        right_value,
                  FROM taxon
                  WHERE taxon_id = %s'''

        name_sql = 'SELECT name, name_class FROM taxon_name WHERE taxon_id = %s'

        taxon_info = self.adaptor.execute_one(taxon_sql, (_id,))
        if not taxon_info:
            raise NotFoundError('No node with id:{0}'.format(_id))
        else:
            name_dict = {}
            name_info = self.adaptor.execute_and_fetchall(name_sql, (_id,))
            for name, name_class in name_info.items():
                try:
                    name_dict[name_class].append(name)
                except KeyError:
                    name_dict[name_class] = [name]

                taxon_name = None
                try:
                    taxon_name = name_dict['scientific name'][0]
                except KeyError:
                    # Maybe warn here that there is no scientific name set
                    pass
            return TaxonNode(self.adaptor, taxon_info[0], taxon_info[2],
                             taxon_info[1], taxon_name, rank=taxon_info[3],
                             genetic_code=taxon_info[4], mito_genetic_code=taxon_info[5],
                             left_val=taxon_info[6], right_val=taxon_info[7], **name_dict)

    def find_elements(self, target=None, terminal=None,
                      ncbi_taxon_id=None, name=None, name_class=None):
        """ Find all nodes in the tree that satisfy the arguments given.

            If there are no arguments given then a ValueError is raised

            :Parameters:
                terminal :
                ncbi_taxon_id : An NCBI taxonomy id
                name : The name of the taxon
                name_class : the type of name. The most common types
                are 'scientific name' or 'synonym'
        """
        if not (ncbi_taxon_id or name or name_class):
            raise ValueError("Please specify at least one of ncbi_taxon_id,"
                             " name or name_class to filter on")

        if name:
            sql = 'SELECT DISTINCT taxon_id FROM taxon_name WHERE name = %s'
            rows = self.adaptor.execute_and_fetch_col0(sql, (name,))
        elif ncbi_taxon_id:
            sql = 'SELECT taxon_id FROM taxon WHERE ncbi_taxon_id = %s'
            rows = self.adaptor.execute_and_fetch_col0(sql, (ncbi_taxon_id,))

        return [self._make_node(x) for x in rows]

    def add(self, name, name_class, rank=None, genetic_code=None,
            mito_genetic_code=None, ncbi_taxon_id=None, parent=None):

        if parent and isinstance(parent, TaxonNode):
            prev = parent._left_val
            self._update_left_right_taxon_values(prev, 2)
        else:
            prev = self.adaptor.execute_one("SELECT MAX(left_value) FROM taxon")[0]

        if not prev:
            prev = 0

        sql = 'INSERT INTO taxon(parent_taxon_id, left_value, right_value,'\
              ' node_rank, mito_genetic_code, genetic_code, ncbi_taxon_id) '\
              'VALUES(%s, %s, %s, %s, %s, %s, %s)'

        self.adaptor.execute(sql, (parent._id, prev + 1, prev + 2, rank,
                                   mito_genetic_code, genetic_code, ncbi_taxon_id))
        taxon_id = self.adaptor.last_id("taxon")

        sql = 'INSERT INTO taxon_name(taxon_id, name, name_class) VALUES(%s, %s, %s)'
        self.adaptor.execute(sql, (taxon_id, name, name_class))

        return self._make_node(taxon_id)

    def remove(self, node):
        if not isinstance(node, TaxonNode):
            raise ValueError("You must pass in a valid TaxonNode object")

        sql = 'DELETE FROM taxon WHERE left_value BETWEEN %s AND %s'
        self.adaptor.execute(sql, (node._left_val, node._right_val))

        self._update_left_right_taxon_values(node._right_val,
                                             node._left_val - node._right_val - 1)

    def move(self, node, parent):
        prev = parent._left_val

        if parent._id == node._id:
            raise IntegrityError('Cannot move node into self')
        if self.is_descendant(node, parent):
            raise IntegrityError('Cannot move node under its own descendant')

        # shift nodes on the width of moving subtree, like in add()
        self._update_left_right_taxon_values(prev, node._right_val - node._left_val + 1)

        sql = 'UPDATE taxon SET parent_taxon_id = %s WHERE taxon_id = %s'
        self.adaptor.execute(sql, (parent._id, node._id))

        # re-fetch node and prev value since could be changed
        node = self._make_node(id)
        target = self.get_node(parent)
        prev = target['left_value']

        sql = '''
            UPDATE taxon
            SET right_value = right_value + %s,
            left_value  = left_value + %s
            WHERE left_value > %s - 1 AND right_value < %s + 1
        '''
        offset = prev - node['left_value'] + 1
        self.adaptor.execute(sql, (offset, offset, node['left_value'], node['right_value']))

        # shift nodes on the width of moving subtree, like in remove()
        self._update_left_right_taxon_values(node['right_value'],
                                             node['left_value'] - node['right_value'] - 1)

    def get_root(self):
        sql = 'SELECT taxon_id FROM taxon WHERE parent_taxon_id IS NULL'\
              ' OR parent_taxon_id = taxon_id'

        results = self.adaptor.execute_and_fetchall(sql)

        if not results:
            raise NotFoundError('No root node')
        elif len(results) > 1:
            raise IntegrityError("Multiple root nodes detected")
        else:
            return self.get_node(results[0][0])


    def get_parent(self, id):
        sql = 'SELECT taxon_id FROM taxon n '\
              'JOIN taxon p ON n.parent_taxon_id = p.taxon_id '\
              'WHERE n.taxon_id = %s'

        result = self.adaptor.execute_one(sql, (id,))

        if not result:
            raise NotFoundError('No parent node for id:{0}'.format(id))
        else:
            return self.get_node(result[0])

    def get_next(self, id):
        sql = 'SELECT taxon_id FROM taxon n '\
              'JOIN taxon s ON n.right_value + 1 = s.left_value '\
              'AND n.parent_taxon_id = s.parent_taxon_id WHERE n.taxon_id = %s'

        result = self.adaptor.execute_one(sql, (id,))

        if not result:
            raise NotFoundError('No right sibling node for id:{0}'.format(id))
        else:
            return self.get_node(result[0])

    def get_previous(self, id):
        sql = 'SELECT taxon_id FROM taxon n '\
              'JOIN taxon s ON n.left_value - 1 = s.right_value '\
              'AND n.parent_taxon_id = s.parent_taxon_id WHERE n.taxon_d = %s'

        result = self.adaptor.execute_one(sql, (id,))

        if not result:
            raise NotFoundError('No left sibling node for id:{0}'.format(id))
        else:
            return self.get_node(result[0])

    def get_path(self, id):
        sql = 'SELECT taxon_id FROM taxon n '\
              'JOIN taxon a ON a.left_value <= n.left_value '\
              'AND a.right_value >= n.right_value '\
              'WHERE n.taxon_id = %s ORDER BY a.left_value ASC'\

        results = self.adaptor.execute_and_fetch_col0(sql, (id,))

        return list(map(self.get_node, results))

    def get_children(self, id):
        sql = 'SELECT taxon_id FROM taxon n '\
              'JOIN taxon c ON n.taxon_id = c.parent_taxon_id '\
              'WHERE n.taxon_id = %s ORDER BY c.left_value ASC'

        results = self.adaptor.execute_and_fetch_col0(sql, (id,))

        return list(map(self.get_node, results))

    def get_descendants(self, id):
        sql = 'SELECT taxon_id FROM taxon d '\
              'JOIN taxon n ON d.left_value BETWEEN n.left_value + 1 '\
              'AND n.right_value - 1 '\
              'WHERE n.taxon_id = %s ORDER BY d.left_value ASC'

        results = self.adaptor.execute_and_fetch_col0(sql, (id,))

        return list(map(self.get_node, results))

    def is_descendant(self, parentId, descendantId):
        parent = self.get_node(parentId)
        descendant = self.get_node(descendantId)

        return parent['left_value'] < descendant['left_value'] and \
               parent['right_value'] > descendant['right_value']

    def isLeaf(self, id):
        node = self.get_node(id)

        return node['right_value'] - node['left_left'] == 1

    def visualize(self, name='name'):
        sql = '''
           SELECT node.taxon_id, node.parent_taxon_id,
           node.left_value, node.right_value,
           taxon_name.name || ' (' || node.node_rank || ')', COUNT(*)
           FROM taxon root
           JOIN taxon node ON node.left_value BETWEEN root.left_value AND root.right_value
           JOIN taxon_name ON node.taxon_id = taxon_name.taxon_id
           WHERE taxon_name.name_class = 'scientific name'
           GROUP BY node.taxon_id, taxon_name.name
           ORDER BY node.left_value ASC
           '''

        results = self.adaptor.execute_and_fetchall(sql)

        return '\n'.join(u'{0}{1} {2} → {3} ({4}, {5})'.format(
            u'  ' * (r[5] - 2) + u'└─' if r[1] else '', r[4], r[0], r[1], r[2], r[3]
            ) for r in results)
